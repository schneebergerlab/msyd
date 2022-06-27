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
import util
import syntools
from syntools import Range
from cigar import Cigar


"""
notes
   - output in file/plot in some clever way
       - maybe as chromosome "bands"?
       - maybe as % syntenic to each query? => clustering/rooted tree
       - use plotsr?

TODO output format

TO/DO handle incorrectly mapped chromosomes (use mapping from syri output)
"""

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering
class Coresyn:
    """
    TODO docstring
    """
    # ranges, cigars have type List[Range]/List[Cigar], respectively, but cython cannot deal with generic type hints
    def __init__(self, ref:Range, ranges, cigars):
        self.ref = ref # optional if using a reference-free algorithm. NONE CURRENTLY IMPLEMENTED!
        self.ranges = ranges
        self.cigars = cigars # length equal to ranges; optional if using approximate matching

    def __repr__(self):
        return f"Coresyn({self.ref}, {self.ranges})"

    def __eq__(l, r):
        return l.ref == r.ref and l.ranges == r.ranges and l.cigars == r.cigars
        
    # for now, only sorts on the reference
    def __lt__(l, r):
        if not l.ref or not r.ref:
            raise ValueError(f"ERROR comparing {l} with {r}: both need to have a reference!")
        return l.ref < r.ref

    def combine(l, r):
        """
        Takes two Coresyn objects and combines them into one, determining the overlap automagically
    .
        At first only implemented to work in a cigar-using, reference-based way.
        TODO implement a cigar-free approximative algorithm
        """

        # compute how much around both sides to drop for each alignment to the reference
        ovstart = max(l.ref.start, r.ref.start)
        ovend = min(l.ref.end, r.ref.end)
        assert(ovstart < ovend, f"ERROR: no overlap found between {l} and {r}")

        rdropstart = ovstart - r.ref.start # 0 if r.ref is maximal, else positive
        ldropstart = ovstart - l.ref.start 

        rdropend = r.ref.end - ovend # 0 if r.ref is minimal, else positive
        ldropend = l.ref.end - ovend

        ref = l.ref.drop(ldropstart, ldropend)
        ranges = [] #TODO switch to dequeue/numpy array
        # calculate the exact position based on the CIGAR strings if both Coresyns have the
        if l.cigars and r.cigars:
            cigars = []

            ## compute the new Ranges and CIGAR strings
            # adjust left Ranges
            for rng, cg in zip(l.ranges, l.cigars):
                start, cg = cg.get_removed(ldropstart, start=True, ref=True)
                end, cg  = cg.get_removed(ldropend, start=False, ref=True)
                ranges.append(rng.drop(start, end))
                cigars.append(cg)

            # adjust right Ranges
            for rng, cg in zip(r.ranges, r.cigars):
                start, cg = cg.get_removed(rdropstart, start=True, ref=True)
                end, cg  = cg.get_removed(rdropend, start=False, ref=True)
                ranges.append(rng.drop(start, end))
                cigars.append(cg)

            return Coresyn(ref, ranges, cigars)
        else:
            # use approximate position calculation
            if l.cigars or r.cigars:
                print(f"WARN: one of {l} or {r} does not have CIGAR strings, falling back to approximate position adjustment!")
                return Coresyn(ref,
                        [rng.drop(ldropstart, ldropend) for rng in l.ranges] +
                        [rng.drop(rdropstart, rdropend) for rng in r.ranges],
                        None)


def find_coresyn(syris, alns, sort=False, ref='a', cores=1):
    """
    Finds core syntenic regions by finding the overlap between all syntenic regions in the input files.
    Fairly conservative.
    :param: a list of filenames of SyRI output and alignment files in BAM, SAM or PAF format to read in, as well as optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be sorted (default False).
    :return: a pandas dataframe containing the chromosome, start and end positions of the core syntenic region for each organism.
    """

    syns = syntools.extract_regions_to_list(syris, ann="SYNAL")

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

    syns = list(map(lambda x: syntools.match_synal(*x, ref=ref), zip(syns, alns)))
    
    # remove overlap
    for syn in syns:
        syntools.remove_overlap(syn)

    #print(syns)

    pansyns = None
    if cores > 1:
        pansyns = util.parallel_reduce(intersect_coresyns, syns, cores)
    else:
        pansyns = functools.reduce(intersect_coresyns, syns)

    return pansyns

def intersect_coresyns(left, right):
    """
    The main business logic of find_coresyn.
    This function takes two dataframes containing syntenic regions and merges each one of them by determining the overlap.
    It also updates the CIGAR string of the alignment accordingly.
    It runs in O(len(left) + len(right)).
    """
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

            if rrow.ref.chr > lrow.ref.chr:
                lrow = next(liter)[1]
                continue
            if lrow.ref.chr > rrow.ref.chr:
                rrow = next(riter)[1]
                continue
            
            # determine if there is an overlap
            ovstart = max(rrow.ref.start, lrow.ref.start)
            ovend = min(rrow.ref.end, lrow.ref.end)
            if ovstart < ovend: # there is valid overlap
                ret.append(Coresyn.combine(lrow, rrow))

            # ratchet by dropping the segment with a smaller end
            if lrow.ref.end > rrow.ref.end: # left is after right
                rrow = next(riter)[1]
            elif rrow.ref.end > lrow.ref.end: # right is after left
                lrow = next(liter)[1]
                # if they stop at the same position, drop the one starting further left
            elif lrow.ref.start > rrow.ref.start:
                rrow = next(riter)[1]
            else: # do whatever
                lrow = next(liter)[1]

        except StopIteration: # nothing more to match
            break

    del riter
    del liter
    return pd.DataFrame(data=ret)




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
#WIP
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
                    # and/or check for intersection
                        byorg[rng.org]= Range(rng.org, rng.haplo, rng.chr, min(rng.start, rngalt.start), max(rng.end, rngalt.end))

                else:
                    byorg[rng.org]= rng
            
            incorporate(qry)
            incorporate(ref)

        ranges.append(pd.DataFrame(byorg, index=[0]))

    return pd.concat(ranges)

def coresyn_from_tsv(path, **kwargs):
    return find_coresyn(*syntools.parse_input_tsv(path), **kwargs)


def collapse_to_df(fins, ref='a', ann="SYN"):
    # Takes a number of input syri files, outputs them all squashed into a single DF
    # For use in graph_pansyn to initialise the graph
    syns = syntools.extract_regions_to_list(fins, ref=ref, ann=ann)
    for syn in syns:
        syn.columns=['ref', 'qry']

    return pd.concat(syns, ignore_index=True)



# make it possible to use this file as test
if __name__ == "__main__": # testing

    import sys

    # removes cigar strings for more concise printing
    remcigar = lambda x: x# x[0] if type(x)==list or type(x)==tuple else x

    syris = []
    alns = []
    for fin in sys.argv[1:]:
        syris.append(fin + "syri.out")
        alns.append(fin + ".bam")

    df1 = find_coresyn(syris, alns, sort=False).apply(lambda x: x.apply(remcigar))
    #print(df1.to_string())
    print("regions:", len(df1))
    print("total lengths:", sum(map(lambda x: x[1][0].end-x[1][0].start, df1.iterrows())))
    sys.exit()
    cg1 = df1.iloc[0, 1][1]
    cg2 = df1.iloc[0, 2][1]
    #print(cg1.to_string())
    #print(cg1.get_identity())
    cgimp = Cigar.impute(cg1, cg2)
    print(cgimp.to_string())
    print(cg1.get_len(), cg1.get_identity(), cg2.get_len(), cg2.get_identity(), cgimp.get_len(), cgimp.get_identity())
