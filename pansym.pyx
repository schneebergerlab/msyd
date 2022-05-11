# -*- coding: utf-8 -*-
#!/usr/bin/python3
# distutils: language = c++
# cython: language_level = 3

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""
import functools
import pandas as pd
import numpy as np
import ingest
#import networkx as nx



# these classes form a part of the general SV format, TODO move them into dedicated file once the format is finalised
# A position is specified by the organism, chromosome, haplotype and base position
# A range takes a start and an end position. If the end < start, the range is defined as inverted
#TODO use proper cython types, e.g. char for haplo
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
    
    def __eq__(l, r):
        return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
                l.start == r.start & l.start == r.start

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


"""
notes
   - TODO add alignment info from bam, use ingest method
        - add CIGAR string extraction & incorporate into position calculation
        - maybe this is the cause of the asymmetry?

   - identify syntenic regions to ref
       - find regions sharing similar locations on the A genome, matching on aStart?
       - sorted join type algorithm?
       - lifting pansyntenic regions? if lifting, what topology?

   - output in file/plot in some clever way
       - maybe as chromosome "bands"?
       - maybe as % syntenic to each query? => clustering/rooted tree
       - use plotsr?

TODO output format

TO/DO be more lenient?
TO/DO handle incorrectly mapped chromosomes (use mapping from syri output)
"""

# refactor out the reading in part, TODO move to ingest and interface to a fileformat-agnostic method once file format finalised
def extract_syn_regions(fins, ref="a"):

    # columns to look for as start/end positions
    reforg = "ref"
    refchr = ref + "chr"
    refhaplo = "NaN"
    refstart = ref + "start"
    refend = ref + "end"

    qry = 'b' if ref == 'a' else 'a' # these seem to be the only two values
    qrychr = qry + "chr"
    qryhaplo = "NaN"
    qrystart = qry + "start"
    qryend = qry + "end"

    syns = []
    for fin in fins:
        qryorg = fin.split('/')[-1]
        buf = []
        raw, chr_mapping = ingest.readsyriout(fin)
        raw = raw.loc[raw['type']=='SYN']
        # if implementing filtering later, filter here

        for row in raw.iterrows():
            row = row[1]
            buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
                Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
                ])

        syns.append(pd.DataFrame(data=buf, columns=["ref", qryorg]))

    return syns


def find_pansyn(fins, ref="a", sort=False):
    """
    Finds pansyntenic regions by finding the overlap between all syntenic regions in the input files.
    Fairly conservative.
    :param: a list of filenames of SyRI output to read in, as well as optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be sorted (default False).
    :return: a pandas dataframe containing the chromosome, start and end positions of the pansyntenic regions on the reference chromosome, as determined by an overlapping intersection
    """

    syns = extract_syn_regions(fins, ref=ref)
    if sort:
        syns = [x.sort_values('ref') for x in syns]

    # take two dataframes of syntenic regions and compute the flexible intersection
    def intersect_syns(left, right):
        #return pd.merge(left, right) # equivalent to old strict intersection
        ret = []

        if len(right) == 0:
            raise ValueError("right is empty!")
        if len(left) == 0:
            raise ValueError("left is empty!")

        riter = right.iterrows()
        rrow = next(riter)[1]
        liter = left.iterrows()
        lrow = next(liter)[1]
        cols = ['ref'] + list(lrow.index)[1:] + list(rrow.index)[1:]
        while True:
            try: # python iterators suck, so this loop is entirely try-catch'ed
                # this may be very slow, TO/DO benchmark # seems to be alright?

                rref = rrow['ref']
                lref = lrow['ref']
                if rref.chr > lref.chr:
                    lrow = next(liter)[1]
                    continue
                if lref.chr > rref.chr:
                    rrow = next(riter)[1]
                    continue
                
                # determine overlap region
                ovstart = max(rref.start, lref.start)
                ovend = min(rref.end, lref.end)

                if ovstart < ovend: # there is valid overlap, i.e. the region is pansyntenic
                    # compute how much around both sides to drop for each alignment to the reference
                    rdropstart = ovstart - rref.start # 0 if rref is maximal, else positive
                    ldropstart = ovstart - lref.start 

                    rdropend = ovend - rref.end # 0 if lref is minimal, else negative
                    ldropend = ovend - lref.end
                    
                    # compute the adjusted position of the pansyntenic regions
                    # drop the references from both columns, place at start (ref is always at first position)
                    ret.append([Range('ref', rref.chr, rref.haplo, ovstart, ovend)] +
                            [Range(syn.org, syn.chr, syn.haplo, syn.start + ldropstart, syn.end + ldropend) for syn in lrow.drop('ref')] +
                            [Range(syn.org, syn.chr, syn.haplo, syn.start + rdropstart, syn.end + rdropend) for syn in rrow.drop('ref')])

                # ratchet by dropping the segment with a smaller end
                if lref.end > rref.end: # left is after right
                    rrow = next(riter)[1]
                elif rref.end > lref.end: # right is after left
                    lrow = next(liter)[1]
                    # if they stop at the same position, drop the one starting further left
                elif lref.start > rref.start:
                    rrow = next(riter)[1]
                else: # do whatever
                    lrow = next(liter)[1]

            except StopIteration: # nothing more to match
                break

        del riter
        del liter
        ret = pd.DataFrame(data=ret, columns=cols)
        return ret.sort_values('ref') if sort else ret


    pansyns = functools.reduce(intersect_syns, syns)

    return pansyns


# use igraph for the graph algorithm backend
# how to encode a region? by organism, coordinates & chromosome?
# how to recognize a region? match over organism, find start position/overlap within tolerances?
# => n^2 algorithm probably
# maybe store nodes sorted? networkx uses a dict by default, but can read from lists
# => do it in a (sorted?) list, then use a similar algorithm to the sequence-based code above to identify syntenic regions?
def exact_pansyn(fins):
    """

    :param:
    :return:
    """
    # helper function to extract only large, syntenic regions suitable for analysis from a DF
    def extract_syn_regions(df):
        return df.loc[df['type']=='SYN']

    syns = [extract_syn_regions(ingest.readsyriout(fin)[0]) for fin in fins]


# make it possible to use this file as test
if __name__ == "__main__": # testing

    import sys

    df = find_pansyn(sys.argv[1:], sort=False)
    print(df)
    print("regions:", len(df))
    print("total lengths:", sum(map(lambda x: x[1]['ref'].end-x[1]['ref'].start,df.iterrows())))
    #df = exact_pansyn(sys.argv[1:])
    #print(df)
    #print("regions:", len(df))
