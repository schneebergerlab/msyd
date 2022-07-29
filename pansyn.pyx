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
import copy

MIN_SYN_THRESH = 0


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

    def __hash__(self):
        return hash(self.org) + hash(self.chr) + hash(self.haplo) + hash(self.pos)

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
        return self.end - self.start + 1 # start is inclusive

    def __hash__(self):
        return hash(self.org) + hash(self.chr) + hash(self.haplo) + hash(self.start) + hash(self.end)

    def get_rightmost(self):
        return max(self.start, self.end)

    def get_leftmost(self):
        return min(self.start, self.end)

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
    The parameter `ranges_dict` is a dictionary of genomic `synctools.Range`es storing the location this syntenic region has on each organism.
    This dictionary cannot be None.
    Other attributes are `ref`, which stores the position on the reference -- this attribute can be `None` if using a reference-free algorithm, but none have been implemented so far.
    `cigars_dict` contains a dictionary of `cigar.Cigar` objects corresponding to the alignment of each `Range` to the reference.
    The attribute can be `None` if using approximate position calculation (usually for performance/memory reasons).\n
    Keys in `cigars_dict` correspond with indices in `ranges_dict`.\n
    In the future, `cigars_dict` may also be used for storing pairwise alignments of the core syntenic regions to each other.
    Also, a separate field for an MSA may be added.

    Pansyn implements comparison operators to enable sorting according to the end on the reference.
    For sorting, the `ref` field needs to be set.

    """
    # ranges_dict, cigars_dict have type Dict[String, Range]/Dict[String, Cigar], respectively, but cython cannot deal with generic type hints
    def __init__(self, ref:Range, ranges_dict, cigars_dict):
        if not ranges_dict:
            raise ValueError(f"ERROR: Trying to initialiase Pansyn with no non-reference Range (ref: {ref})")
        if cigars_dict and not ranges_dict.keys() == cigars_dict.keys():
            raise ValueError(f"ERROR: Trying to initialise Pansyn with ranges_dict keys {ranges_dict.keys()} not matching cigars_dict keys {cigars_dict.keys()}!")
        self.ref = ref # optional if using a reference-free algorithm. NONE CURRENTLY IMPLEMENTED!
        self.ranges_dict = ranges_dict
        self.cigars_dict = cigars_dict # optional if using approximate matching

    def __repr__(self):
        return f"Pansyn({self.ref}, {self.ranges_dict})"

    def __eq__(l, r):
        return l.ref == r.ref and l.ranges_dict == r.ranges_dict and l.cigars_dict == r.cigars_dict
        
    # for now, only sorts on the reference (falling back to the Range comparison operator)
    def __lt__(l, r):
        if not l.ref or not r.ref:
            raise ValueError(f"ERROR comparing {l} with {r}: both need to have a reference!")
        return l.ref < r.ref

    def __hash__(self):
        return hash(self.ref)# + hash(self.ranges_dict) + hash(self.cigars_dict) # caused problems with deque

    def add(self, rng:Range, cg: Cigar):
        self.ranges_dict[rng.org] = rng
        if cg:
            if self.cigars_dict:
                self.cigars_dict[rng.org] = cg
            else:
                print("WARNING: attempted to add cigar to Pansyn without cigars_dict, ignoring")

    def degree(self):
        return len(self.ranges_dict)

    def get_organisms(self):
        return self.ranges_dict.keys()

    def get_ranges(self):
        return self.ranges_dict.values()

    def __add__(self, other):
        """
        Convenience function to concatenate two `Pansyn` objects.
        Uses a shallow copy of the cigar/range to stay without side effects.
        """
        rngs = copy.copy(self.ranges_dict)
        rngs.update(other.ranges_dict)

        cgs = None
        if self.cigars_dict and other.cigars_dict:
            cgs = copy.copy(self.cigars_dict)
            cgs.update(other.cigars_dict)
        elif self.cigars_dict or other.cigars_dict:
            print(f"WARN: Trying to add two Pansyns {self}, {other} with one having CIGARs and one not! Discarding CIGARS!")

        return Pansyn(self.ref, rngs, cgs)

    def drop(self, start, end):
        """
        Returns a new `Pansyn` object with `start`/`end` positions from the start/end of this pansyntenic region removed, respecting cigar alignments if not `None`.
        """
        oldref = self.ref
        ref = self.ref.drop(start, end)
        ranges_dict = dict()
        cigars_dict = None
        if not self.cigars_dict:
            ranges_dict = {org:rng.drop(start, end) for (org, rng) in self.ranges_dict.items()}
        else:
            cigars_dict = dict()
            for org, rng in self.ranges_dict.items():
                cg = self.cigars_dict[org]
                try:
                    start, cg = cg.get_removed(start, start=True, ref=True)
                    end, cg = cg.get_removed(end, start=False, ref=True)
                except ValueError:
                    print(f"ERROR: invalid input to cg.get_removed({start}, {end}) on {rng} (len: {len(rng)}). Check if start, end are correct!")
                    if len(rng) < 2000:
                        print("ref is:", ref, "was:", oldref)
                        print("cigar was:", self.cigars_dict[org])
                    continue
                ranges_dict[org] = rng.drop(start, end)
                cigars_dict[org] = cg

        return Pansyn(ref, ranges_dict, cigars_dict)



def calc_overlap(l: Pansyn, r: Pansyn, detect_crosssyn=False, allow_overlap=False):
    """
    Used to be called `combine` in coresyn and crosssyn.
    Takes two `Pansyn` objects, finds the overlap and returns a list of the `Pansyn` objects generated by overlapping the two.
    Returns a sorted list of at most three `Pansyn` objects.
    :param: `detect_crosssyn` specifies whether to output anything but the overlap.
    If disabled will only find core syntenic regions.
    `allow_overlap` specifies whether pansyntenic regions are allowed to overlap, in essence outputting the input regions as-is.
    Output may be non-associative and non-commutative if set to `True`.
    """

    if not l.ref and not r.ref:
        raise ValueError("ERROR: calc_overlap can only be called on two Pansyns with references!")

    # compute overlap
    ovstart = max(l.ref.start, r.ref.start)
    ovend = min(l.ref.end, r.ref.end)
    if not ovstart < ovend:
        print(f"WARNING: no overlap found between {l} and {r}, returning empty!")
        return []

    leftest = l if l.ref.start < r.ref.start else r
    rightest = l if l.ref.end > r.ref.end else r

    ret = set() # use a set to automatically remove duplicates

    def add_lenfiltered(pansyn): # checks if the region is higher than MIN_SYN_THRESH, if so, adds it to ret
        if len(pansyn.ref) > MIN_SYN_THRESH:
            ret.add(pansyn)

    if detect_crosssyn:
        if allow_overlap:
            add_lenfiltered(leftest)
        else:
            add_lenfiltered(leftest.drop(0, leftest.ref.end - ovstart))
    
    # core synteny
    add_lenfiltered(l.drop(ovstart - l.ref.start, l.ref.end - ovend) + r.drop(ovstart - r.ref.start, r.ref.end - ovend))
    # this sometimes tries to drop more than a range has – maybe handle this properly somehow?
    # idea from manish: check if cigar string length also corresponds to length in query, to check for errors in code/data

    if detect_crosssyn:
        if allow_overlap:
            add_lenfiltered(rightest)
        else:
            add_lenfiltered(rightest.drop(rightest.ref.start - ovend, 0))

    return sorted(ret)


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
            org = synr[1].org
            if synr[0].chr == alnr[refchr] and synr[0].start == alnr[refstart] and synr[0].end == alnr[refend]:
                ret.append(Pansyn(ref=synr[0], ranges_dict={org:synr[1]}, cigars_dict={org:Cigar.from_string(alnr['cg'])}))
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
        if cur.cigars_dict:
            for org, cg in cur.cigars_dict.items():
                drop, cgnew = cg.get_removed(ov, start=True)
                cur.cigars_dict[org] = cgnew
                cur.ranges_dict[org].start += drop
        else:
            for rng in cur.ranges_dict.values():
                rng.start += ov

        prev = cur


def find_overlaps(left, right, **kwargs):
    """
    This function takes two dataframes containing syntenic regions and outputs the overlap found between each of them as a new pandas dataframe.
    It runs in O(len(left) + len(right)).
    Calls `calc_overlap` to determine the actual overlaps to merge as output.
    """
    ret = deque()

    if len(right) == 0:
        raise ValueError("right is empty!")
    if len(left) == 0:
        raise ValueError("left is empty!")

    riter = right.iterrows()
    rrow = next(riter)[1][0]
    liter = left.iterrows()
    lrow = next(liter)[1][0]
    while True:
        try: # python iterators suck, so this loop is entirely try-catch'ed
            lref = lrow.ref
            rref = rrow.ref
            if rref.chr > lref.chr:
                lrow = next(liter)[1][0]
                continue
            if lref.chr > rref.chr:
                rrow = next(riter)[1][0]
                continue
            
            # determine if there is an overlap
            ovstart = max(rref.start, lref.start)
            ovend = min(rref.end, lref.end)
            if ovend - ovstart > MIN_SYN_THRESH: # there is valid overlap
                ret.extend(calc_overlap(lrow, rrow, **kwargs))

            # ratchet by dropping the segment with a smaller end
            if lref.end > rref.end: # left is after right
                rrow = next(riter)[1][0]
            elif rref.end > lref.end: # right is after left
                lrow = next(liter)[1][0]
                # if they stop at the same position, drop the one starting further left
            elif lref.start > rref.start:
                rrow = next(riter)[1][0]
            else: # do whatever
                lrow = next(liter)[1][0]

        except StopIteration: # nothing more to match
            break

    del riter
    del liter
    return pd.DataFrame(data=list(ret))


def find_multisyn(syris, alns, sort=False, ref='a', cores=1, SYNAL=False, **kwargs):
    """
    Finds core and cross-syntenic regions in the input files, depending on if the parameter `detect_crossyn` that is ultimately passed on to `calc_overlap` is set to `True`.
    Fairly conservative.
    :param: a list of filenames of SyRI output and alignment files in BAM, SAM or PAF format to read in, parameters optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be sorted (default False).
    `alns` can be set to `None`, in which case the region lengts will be estimated instead of calculated exactly from CIGAR strings.
    Crosssynteny is detected if the `detect_crosssyn` parameter is set to `True`.
    :return: a pandas dataframe containing the chromosome, start and end positions of the core syntenic region for each organism.
    """

    syns = extract_regions_to_list(syris, ann="SYNAL" if SYNAL else "SYN")

    if sort:
        syns = [x.sort_values(x.columns[0]) for x in syns]

    alnfilelookup = {
            'sam': ingest.readSAMBAM,
            'bam': ingest.readSAMBAM,
            'paf': ingest.readPAF
            }

    if alns and not SYNAL:
        alns = [alnfilelookup[aln.split('.')[-1]](aln) for aln in alns]
        alns = [aln[(aln.adir==1) & (aln.bdir==1)] for aln in alns] # only count non-inverted alignments as syntenic
        #print(alns)

        syns = [match_synal(*x, ref=ref) for x in zip(syns, alns)]
    else:
        syns = [
                pd.DataFrame(
                    [ Pansyn(ref=row[1][0],
                        ranges_dict={row[1][1].org:row[1][1]})
                        for row in s.iterrows()]) for s in syns]

    #for syn in syns:
    #    for pansyn in syn.iterrows():
    #        pansyn = pansyn[1][0]
    #        cg = list(pansyn.cigars_dict.values())[0]
    #        rng = list(pansyn.ranges_dict.values())[0]
    #        # at this point, each pansyn should only have a reference and one query region
    #        if len(pansyn.ref) != cg.get_len(ref=True):
    #            print(f"ERRORERRROR: cigar string length {cg.get_len(ref=True)} does not match reference length in {pansyn}")
    #            #print(cg)
    #        if len(rng) != cg.get_len(ref=False):
    #            print(f"ERRORERRROR: cigar string length {cg.get_len(ref=False)} does not match query length in {pansyn}")
    #            print(cg)
    #TODO weird bug: WTF, the lengths of cigar/alignment do not match???
    #TODO weird bug: maybe mistake in reading in from alignment file?? maybe double-check alignments for correctness

    # remove overlap
    for syn in syns:
        remove_overlap(syn)
    print("INFO: overlap removed")

    pansyns = None
    ovlap = functools.partial(find_overlaps, **kwargs)
    if cores > 1:
        pansyns = util.parallel_reduce(ovlap, syns, cores)
    else:
        pansyns = functools.reduce(ovlap, syns)

    return pansyns

