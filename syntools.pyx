#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import numpy as np
import ingest
import util
from cigar import Cigar
import functools
from collections import deque


# these classes form a part of the general SV format
# A position is specified by the organism, chromosome, haplotype and base position
# A range takes a start and an end position. If the end < start, the range is defined as inverted
#TODO use proper cython types, e.g. char for haplo

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
class Position:
    def __init__(self, org:str, chr:int, haplo:str, pos: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.pos = pos

    def __repr__(self):
        return f"Position({self.org}, {self.chr}, {self.haplo}, {self.pos})"

    def __eq__(l, r):
        return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
                l.pos == r.pos

    def __lt__(l, r):
        if l.org != r.org:
            raise ValueError("Comparison between different organisms!")
        if l.chr < r.chr:
            return True
        elif l.chr == r.chr:
            return l.pos < r.pos
        else:
            return False

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
class Range:
    def __init__(self, org:str, chr:int, haplo:str, start: int, end: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.start = start
        self.end = end

    def __repr__(self):
        return f"Range({self.org}, {self.chr}, {self.haplo}, {self.start}, {self.end})"
    
    #def __eq__(l, r):
    #    return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
    #            l.start == r.start & l.start == r.start

    # this operator sorts according to the END, not start value,
    # to enable the end ratchet to work properly
    # TO/DO possible refactor: sort by start here, invert in sorting for algorithm
    # shouldn't really matter as the regions are nonoverlapping, but...
    def __lt__(l, r):
        if l.org != r.org:
            raise ValueError("Comparison between different organisms!")

        if l.chr < r.chr:
            return True
        elif l.chr == r.chr:
            if l.end < r.end:
                return True
            elif l.end == r.end:
                return l.start < r.start
        return False

    def __len__(self):
        return self.start - self.end

    def drop(self, start, end):
        """
        :param: 'start'/'end' specify how much to drop on each end.
        :return: A Range missing the specified area.
        """
        return Range(self.org, self.chr, self.haplo, self.start + start, self.end - end)


# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering
class Pansyn:
    """
    A class representing a region syntenic among a set of genomes.
    The parameter `ranges` is a list of genomic `synctools.Range`es storing the locations this syntenic regions has on each organism.
    This list cannot be None and should contain an entry for every genome being compared.
    Other attributes are `ref`, which stores the position on the reference -- this attribute can be `None` if using a reference-free algorithm, but none have been implemented so far.
    `cigars` contains a list of `cigar.Cigar` objects corresponding to the alignment of each position to the reference.
    The attribute can be `None` if using approximate position calculation (usually for performance/memory reasons).\n
    Indices in `cigars` correspond with indices in `ranges`.
    In the future, `cigars` may also be used for storing pairwise alignments of the core syntenic regions to each other.
    Also, a separate field for an MSA may be added.

    Pansyn implements comparison operators to enable sorting according to the end on the reference.
    For sorting, the `ref` field needs to be set.

    """
    # ranges, cigars have type List[Range]/List[Cigar], respectively, but cython cannot deal with generic type hints
    def __init__(self, ref:Range, ranges, cigars):
        self.ref = ref # optional if using a reference-free algorithm. NONE CURRENTLY IMPLEMENTED!
        self.ranges = ranges
        self.cigars = cigars # length equal to ranges; optional if using approximate matching

    def __repr__(self):
        return f"Pansyn({self.ref}, {self.ranges})"

    def __eq__(l, r):
        return l.ref == r.ref and l.ranges == r.ranges and l.cigars == r.cigars
        
    # for now, only sorts on the reference (falling back to the Range comparison operator)
    def __lt__(l, r):
        if not l.ref or not r.ref:
            raise ValueError(f"ERROR comparing {l} with {r}: both need to have a reference!")
        return l.ref < r.ref

    def add(self, rng:Range, cg: Cigar):
        self.ranges.append(rng)
        self.cigars.append(cg)

    def degree(self):
        return len(self.ranges)

    def __add__(self, other):
        """
        Convenience function to concatenate two `Pansyn` objects
        """
        return Pansyn(self.ref, self.ranges + other.ranges, self.cigars + other.cigars)

    def drop(self, start, end):
        """
        Rremoves `start`/`end` positions from the start/end of this pansyntenic region, respecting cigar alignments if not `None`.
        """
        ret_rngs = deque()
        ret_cgs = deque()
        if not self.cigars:
            self.ref = self.ref.drop(start, end)
            self.ranges = [rng.drop(start, end) for rng in self.ranges]
            return
        for i in range(len(self.ranges)):
            start, cg = self.cigars[i].get_removed(start, start=True, ref=True)
            end, cg = cg.get_removed(end, start=False, ref=True)
            self.ranges[i] = self.ranges[i].drop(start, end)
            self.cigars[i] = cg

        return (ret_rngs, ret_cgs)


# given a bam file and corresponding SYNAL range df,
# Transform them into one list of Pansyn objects
def match_synal(syn, aln, ref='a'):
    """
    This function takes an aligment and SYNAL dataframe and matches corresponding regions.
    It returns a dataframe containing the regions with the corresponding CIGAR string as a `Pansyn` object.
    :params: syn: SYNAL dataframe, aln: alignment dataframe, ref: whether the reference is the 'a' or 'b' strand in the alignment dataframe.
    :returns: a dataframe containing the SYNAL regions with corresponding CIGAR strings as `Pansyn` objects.
    """
    ret = deque()
    syniter = syn.iterrows()
    alniter = aln.iterrows()
    refchr = ref + "chr"
    refstart = ref + "start"
    refend = ref + "end"

    synr = next(syniter)[1]
    alnr = next(alniter)[1]
    while True:
        try:
            if synr[0].chr == alnr[refchr] and synr[0].start == alnr[refstart] and synr[0].end == alnr[refend]:
                ret.append(Pansyn(ref=synr[0], ranges=[synr[1]], cigars=[Cigar.from_string(alnr['cg'])]))
                synr = next(syniter)[1]
            alnr = next(alniter)[1]
        except StopIteration:
            break

    return pd.DataFrame(list(ret))

def extract_regions(fin, ref='a', ann='SYN', reforg='ref', qryorg='qry'):
    """
    Given a syri output file, extract all regions matching a given annotation.
    """

    # columns to look for as start/end positions
    refchr = ref + "chr"
    refhaplo = "NaN"
    refstart = ref + "start"
    refend = ref + "end"

    qry = 'b' if ref == 'a' else 'a' # these seem to be the only two values in syri output
    qrychr = qry + "chr"
    qryhaplo = "NaN"
    qrystart = qry + "start"
    qryend = qry + "end"

    buf = deque()
    raw, chr_mapping = ingest.readsyriout(fin) #TODO? handle chr_mapping
    raw = raw.loc[raw['type'] == ann]
    # if implementing filtering later, filter here

    for row in raw.iterrows():
        row = row[1]
        buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
            Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
            ])

    return pd.DataFrame(data=list(buf), columns=[reforg, qryorg])

def extract_regions_to_list(fins, **kwargs):
    """
    `extract_regions`, but for processing a list of inputs
    """
    return [extract_regions(fin, **kwargs,\
            reforg=fin.split('/')[-1].split('_')[0],\
            qryorg=fin.split('/')[-1].split('_')[-1].split('syri')[0])\
            for fin in fins]


def remove_overlap(syn):
    """
    part of the preprocessing of SYNAL regions for find_multisyn
    removes overlap from the first region if two overlapping regions are next to each other
    assumes syn to be sorted
    mutates syn
    """
    syniter = syn.iterrows()
    prev = next(syniter)[1][0]
    for _, cur in syniter:
        cur = cur[0]
        # check for overlap on the reference
        ov = prev.ref.end - cur.ref.start
        if ov <= 0 or cur.ref.chr != prev.ref.chr: # there is no overlap
            prev = cur
            continue

        #print(ov, prevref, curref)
        # remove the overlap from the previous row;
        # this performs essentially the same function as compute_overlap
        # in intersect_syns
        cur.ref.start += ov
        for ind, cg in enumerate(cur.cigars):
            drop, cgnew = cg.get_removed(ov, start=True)
            cur.cigars[ind] = cgnew
            cur.ranges[ind].start += drop

        prev = cur


def parse_input_tsv(path):
    """
    Takes a file containing the input alignments/syri files and processes it for find_multisyn.
    Anything after a # is ignored. Lines starting with # are skipped.
    :params: path to a file containing the paths of the input alignment and syri files in tsv format
    :returns: a tuple of two lists containing the paths of the alignment and syri files.
    """
    import os
    syris = deque()     # Lists are too slow appending, using deque instead
    alns = deque()
    with open(path, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue

            val = line.strip().split('#')[0].split('\t')
            if len(val) > 2:
                print(f"ERROR: invalid entry in {path}. Skipping line: {line}")
                continue
            # Check that the files are accessible
            if not os.path.isfile(val[0]):
                raise FileNotFoundError(f"Cannot find file at {val[0]}. Exiting")
            if not os.path.isfile(val[1]):
                raise FileNotFoundError(f"Cannot find file at {val[1]}. Exiting")

            alns.append(val[0].strip())
            syris.append(val[1].strip())

    return (syris, alns)


def find_multisyn(syris, alns, intersect, sort=False, ref='a', cores=1):
    """
    Finds core or cross-syntenic regions in the input files, depending on the intersection operation specified in `intersect`.
    Fairly conservative.
    :param: a list of filenames of SyRI output and alignment files in BAM, SAM or PAF format to read in, the intersection operation and class to use (from either `coresyn` or `crossyn`) as well as optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be sorted (default False).
    :return: a pandas dataframe containing the chromosome, start and end positions of the core syntenic region for each organism.
    """

    syns = extract_regions_to_list(syris, ann="SYNAL")

    if sort:
        syns = [x.sort_values(x.columns[0]) for x in syns]

    alnfilelookup = {
            'sam': ingest.readSAMBAM,
            'bam': ingest.readSAMBAM,
            'paf': ingest.readPAF
            }

    alns = [alnfilelookup[aln.split('.')[-1]](aln) for aln in alns]
    alns = [aln[(aln.adir==1) & (aln.bdir==1)] for aln in alns] # only count non-inverted alignments as syntenic
    #print(alns)

    syns = list(map(lambda x: match_synal(*x, ref=ref), zip(syns, alns)))
    
    # remove overlap
    for syn in syns:
        remove_overlap(syn)

    #print(syns)

    pansyns = None
    #TODO switch to crosssyn
    if cores > 1:
        pansyns = util.parallel_reduce(intersect, syns, cores)
    else:
        pansyns = functools.reduce(intersect, syns)

    return pansyns

