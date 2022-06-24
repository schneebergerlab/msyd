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
from syn import Range
from cigar import Cigar


"""
notes
   - output in file/plot in some clever way
       - maybe as chromosome "bands"?
       - maybe as % syntenic to each query? => clustering/rooted tree
       - use plotsr?

TODO output format

TO/DO be more lenient? Store "support" for each region as pct or number having/not having this?
=> pleio/polysynteny
TO/DO handle incorrectly mapped chromosomes (use mapping from syri output)
"""


def collapse_to_df(fins, ref='a', ann="SYN"):
    # Takes a number of input syri files, outputs them all squashed into a single DF
    # For use in graph_pansyn to initialise the graph
    syns = syntools.extract_regions_to_list(fins, ref=ref, ann=ann)
    for syn in syns:
        syn.columns=['ref', 'qry']

    return pd.concat(syns, ignore_index=True)

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
        pansyns = util.parallel_reduce(intersect_syns, syns, cores)
    else:
        pansyns = functools.reduce(intersect_syns, syns)

    return pansyns

def intersect_syns(left, right):
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

                rdropend = rref.end - ovend # 0 if rref is minimal, else positive
                ldropend = lref.end - ovend
                
                # function for computing the adjusted position of the pansyntenic regions
                def compute_position(tup, dropstart, dropend):
                    syn = tup[0]
                    cg = tup[1]

                    #print(syn, dropstart, dropend)
                    #from collections import defaultdict
                    #stats = defaultdict(lambda: 0)
                    #for x in cg:
                    #    stats[x[1]] = stats[x[1]] + x[0]
                    #print(stats)

                    start, cg = cg.get_removed(dropstart, start=True, ref=True)
                    end, cg  = cg.get_removed(dropend, start=False, ref=True)
                    return [syn.drop(start, end), cg]

                # drop the references from both rows, place at start (ref is always at first position)
                ret.append([rref.drop(rdropstart, rdropend)] +
                        [compute_position(tup, ldropstart, ldropend) for tup in lrow.drop(lrow.index[0])] +
                        [compute_position(tup, rdropstart, rdropend) for tup in rrow.drop(rrow.index[0])])

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
    return pd.DataFrame(data=ret, columns=cols)
    #return ret.sort_values(ret.columns[0]) if sort else ret




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

def parse_input_tsv(path):
    """
    Takes a file containing the input alignments/syri files and processes it for coresyn.pyx.
    Anything after a # is ignored. Lines starting with # are skipped.
    :params: path to a file containing the paths of the input alignment and syri files in tsv format
    :returns: a tuple of two lists containing the paths of the alignment and syri files.
    """
    from collections import deque
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

def coresyn_from_tsv(path, **kwargs):
    return find_coresyn(*parse_input_tsv(path), **kwargs)



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
