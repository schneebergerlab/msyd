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
import networkx as nx


# notes
#   - start with syri output, add alignment info later if needed
#   - identify syntenic regions to ref
#       - find regions sharing similar locations on the A genome, matching on aStart?
#       - sorted join type algorithm?
#       - lifting pansyntenic regions? if lifting, what topology?
#   - output in file/plot in some clever way
#       - maybe as chromosome "bands"?
#       - maybe as % syntenic to each query? => clustering/rooted tree
#       - use plotsr?

def find_pansyn(fins, ref="a"):
    """
    Finds pansyntenic regions by finding the overlap between all syntenic regions in the input files.
    Seems to be very conservative.
    :param: a list of filenames of SyRI output to read in, as well as optionally specifying which sequence is the reference (default 'a').
    :return: a pandas dataframe containing the chromosome, start and end positions of the pansyntenic regions on the reference chromosome, as determined by an overlapping intersection
    """
    # columns to look for as start/end positions
    refchr = ref + "chr"
    refstart = ref + "start"
    refend = ref + "end"

    # helper function to extract only large, syntenic regions suitable for analysis from a DF
    def extract_syn_regions(df):
        return df.loc[df['type']=='SYN'][[refchr, refstart, refend]]

    syns = [extract_syn_regions(ingest.readsyriout(fin)[0]) for fin in fins]

    # take two dataframes of syntenic regions and compute the flexible intersection
    def intersect_syns(left, right):
        #return pd.merge(left, right) # equivalent to old strict intersection
        
        # they should already be sorted -- required for the rest of the algorithm
        #left.sort_values([refchr, refstart, refend])
        #right.sort_values([refchr, refstart, refend])

        ret = pd.DataFrame(columns = left.columns) # or delete from left? might be faster
        #ret = []

        if len(right) == 0:
            raise ValueError("right is empty!")
        if len(left) == 0:
            raise ValueError("left is empty!")

        riter = right.iterrows()
        rsyn = next(riter)[1]
        liter = left.iterrows()
        lsyn = next(liter)[1]
        ldropped = 0
        rdropped = 0
        while True:
            try: # python iterators suck, so this loop is entirely try-catch'ed
                # this may be very slow, TODO benchmark

                # this uses lexicalic comparisons on strings and is most likely fairly slow
                # TODO possible improvement: store chromosomes as unsigned byte (cython)
                if rsyn[refchr] > lsyn[refchr]:
                    lsyn = next(liter)[1]
                    ldropped = ldropped + 1
                    continue
                if lsyn[refchr] > rsyn[refchr]:
                    rsyn = next(riter)[1]
                    rdropped = rdropped + 1
                    continue
                
                # determine overlap region
                ovstart = max(rsyn[refstart], lsyn[refstart])
                ovend = min(rsyn[refend], lsyn[refend])

                if ovstart < ovend: # there is valid overlap
                    ret = pd.concat([ret, pd.DataFrame(data=[[lsyn[refchr], ovstart, ovend]], columns=left.columns)])
                    # if an overlap has been found, the segments will not be dropped
                    rdropped = rdropped -1
                    ldropped = ldropped -1

                # ratchet by dropping the segment with a smaller start
                if lsyn[refstart] > rsyn[refstart]: # left is after right
                    rsyn = next(riter)[1]
                    rdropped = rdropped + 1
                else: # right is after left
                    lsyn = next(liter)[1]
                    ldropped = ldropped + 1
            except StopIteration: # nothing more to match
                break

        del riter
        del liter
        print("dropped l:", ldropped, "/", len(left), "r:", rdropped, "/", len(right))
        return ret


    pansyns = functools.reduce(intersect_syns, syns)

    return pansyns


# use networkx for the graph algorithm backend
# how to encode a region? by organism, coordinates & chromosome?
#  
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
    df = find_pansyn(sys.argv[1:])
    print(df)
    print("regions:", len(df))
    print("total lengths:", sum(map(lambda x: x[1][2]-x[1][1],df.iterrows())))
