#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""
import functools
import pandas as pd
import numpy as np
import ingest



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
def extract_regions(fin, ref='a', ann='SYN', reforg='ref', qryorg='qry'):

    # columns to look for as start/end positions
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
    buf = []
    raw, chr_mapping = ingest.readsyriout(fin)
    raw = raw.loc[raw['type'] == ann]
    # if implementing filtering later, filter here

    for row in raw.iterrows():
        row = row[1]
        buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
            Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
            ])

    return pd.DataFrame(data=buf, columns=[reforg, qryorg])

def extract_regions_to_list(fins, ref='a', ann="SYN"):
    return [extract_regions(fin, ann=ann,\
            reforg=fin.split('/')[-1].split('_')[0],\
            qryorg=fin.split('/')[-1].split('_')[-1].split('syri')[0])\
            for fin in fins]

def collapse_to_df(fins, ref='a', ann="SYN"):
    syns = extract_regions_to_list(fins, ref=ref, ann=ann)
    for syn in syns:
        syn.columns=['ref', 'qry']

    return pd.concat(syns, ignore_index=True)


def find_pansyn(fins, sort=False):
    """
    Finds pansyntenic regions by finding the overlap between all syntenic regions in the input files.
    Fairly conservative.
    :param: a list of filenames of SyRI output to read in, as well as optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be sorted (default False).
    :return: a pandas dataframe containing the chromosome, start and end positions of the pansyntenic region for each organism.
    """

    syns = extract_regions_to_list(fins)
    print(syns)
    if sort:
        syns = [x.sort_values(x.columns[0]) for x in syns]

    # take two dataframes of syntenic regions and compute the flexible intersection
    def intersect_syns(left, right):
        ret = []

        if len(right) == 0:
            raise ValueError("right is empty!")
        if len(left) == 0:
            raise ValueError("left is empty!")

        riter = right.iterrows()
        rrow = next(riter)[1]
        liter = left.iterrows()
        lrow = next(liter)[1]
        cols = list(lrow.index) + list(rrow.index)[1:]
        while True:
            try: # python iterators suck, so this loop is entirely try-catch'ed
                # this may be very slow, TO/DO benchmark # seems to be alright?

                rref = rrow[0] # the reference always is stored first
                lref = lrow[0]
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
                    ret.append([Range(rref.org, rref.chr, rref.haplo, ovstart, ovend)] +
                            [Range(syn.org, syn.chr, syn.haplo, syn.start + ldropstart, syn.end + ldropend) for syn in lrow.drop(lrow.index[0])] +
                            [Range(syn.org, syn.chr, syn.haplo, syn.start + rdropstart, syn.end + rdropend) for syn in rrow.drop(rrow.index[0])])

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
        return ret.sort_values(ret.columns[0]) if sort else ret


    pansyns = functools.reduce(intersect_syns, syns)

    return pansyns



"""
- use igraph for the graph algorithm backend
- adapt the overlap-algorithm:
    - insert sorted Ranges into igraph, Range object as decorator
    - different levels of strictness: do an edge between two nodes if:
        - subregion: one contains the other
        - overlap: one has an overlap with the other
        - tolerance: their start/end are within a certain tolerance
- find cliques, if that doesn't work try clusters
- then extract information from cliques
"""
#TODO use pairwise syri information instead of reference, match on organism
#TODO try to get sorted algorithm to work properly
def graph_pansyn(fins, mode="overlap", tolerance=100):
    """

    :param:
    :return:
    """
    import igraph

    ## prepare & load the vertices
    colsyns = collapse_to_df(fins)
    colsyns = colsyns.sort_values('ref')

    g = igraph.Graph()
    g.add_vertices(len(colsyns))
    g.vs["ref"] = colsyns['ref']
    g.vs["qry"] = colsyns['qry']

    # try brute force, the algo below seems to have some problems
    for a in g.vs:
        for b in g.vs:
            aref = a['ref']
            bref = b['ref']
            aqry = a['qry']
            bqry = b['qry']

            if not aref.chr == bref.chr:
                continue
            if mode == 'overlap':
                if min(aref.end, bref.end) >= max(aref.start, bref.end):
                    g.add_edge(a.index, b.index)
                if min(aqry.end, bqry.end) >= max(aqry.start, bqry.end):
                    g.add_edge(a.index, b.index)


            elif mode == 'tolerance':
                if abs(aref.end - bref.end) < tolerance and abs(aref.start - bref.start) < tolerance:
                    g.add_edge(a.index, b.index)
                if abs(aqry.end - bqry.end) < tolerance and abs(aqry.start - bqry.start) < tolerance:
                    g.add_edge(a.index, b.index)
                
            elif mode == 'subregion':
                if (aref.start < bref.start) == (aref.end > bref.end):
                    g.add_edge(a.index, b.index)
                if (aqry.start < bqry.start) == (aqry.end > bqry.end):
                    g.add_edge(a.index, b.index)
            


    ## add the appropriate edges, applying the criterion specifyed in mode
#    aiter = iter(g.vs)
#    a = next(aiter)
#    biter = iter(g.vs)
#    b = next(biter)
#    while True:
#        try:
#            aref = a['ref']
#            bref = b['ref']
#            if aref.chr > bref.chr:
#                a = next(aiter)
#                continue
#            elif aref.chr < bref.chr:
#                b = next(biter)
#                continue
#
#            if mode == 'overlap':
#                if min(aref.end, bref.end) >= max(aref.start, bref.end):
#                    g.add_edge(a.index, b.index)
#
#            elif mode == 'tolerance':
#                if abs(aref.end - bref.end) < tolerance and abs(aref.start - bref.start) < tolerance:
#                    g.add_edge(a.index, b.index)
#                
#            elif mode == 'subregion':
#                if (aref.start < bref.start) == (aref.end > bref.end):
#                    g.add_edge(a.index, b.index)
#                
#            else:
#                raise ValueError(f"Invalid mode {mode}")
#            
#            # ratchet by dropping the segment with a smaller end
#            if aref.end > bref.end: # left is after right
#                b = next(biter)
#            elif bref.end > aref.end: # right is after left
#                a = next(aiter)
#            # if they stop at the same position, drop the one starting further left
#            elif aref.start > bref.start:
#                b = next(biter)
#            else: # do whatever
#                a = next(aiter)
#
#        except StopIteration: # nothing to match
#            break
#    del aiter
#    del biter

    cliques = g.maximal_cliques()
    #clusters = g.community_leiden()
    
    ## extract the syntenic cliques, transform them to output

    ranges = []
    for cl in cliques:
        byorg = {}
        for node in cl:
            qry = g.vs[node]['qry']
            ref = g.vs[node]['ref']
            # helper function to incorporate node ranges into the dict
            def incorporate(rng):
                if rng.org in byorg:
                    rngalt = byorg[rng.org]

                    
                    if rngalt.chr == rng.chr:
                        #raise ValueError(f"Mismatching chromosomes unioning {rng} & {rngalt}")
                    # resolve multiple ranges on the same organism as the union of the two
                    # is this appropriate? maybe use intersection
                    # and/oder check for intersection
                        byorg[rng.org]= Range(rng.org, rng.haplo, rng.chr, min(rng.start, rngalt.start), max(rng.end, rngalt.end))

                else:
                    byorg[rng.org]= rng
            
            incorporate(qry)
            incorporate(ref)

        ranges.append(pd.DataFrame(byorg, index=[0]))

    return pd.concat(ranges)




# make it possible to use this file as test
if __name__ == "__main__": # testing

    import sys

    df = find_pansyn(sys.argv[1:], sort=False)
    print(df)
    print("regions:", len(df))
    print("total lengths:", sum(map(lambda x: x[1][0].end-x[1][0].start,df.iterrows())))
    df = graph_pansyn(sys.argv[1:], mode='overlap')
    print(df)
    print("regions:", len(df))
