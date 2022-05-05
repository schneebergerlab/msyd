# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
#!/usr/bin/python3

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""
import functools
import pandas as pd
import numpy as np
import ingest


# notes
#   - start with syri output, add alignment info later if needed
#   - identify syntenic regions to ref, possibly reusing code from syri(synsearch.pyx, ll. 500-550)
#       - find regions sharing similar locations on the A genome, matching on aStart?
#       - sorted join type algorithm?
#       - lifting pansyntenic regions? if lifting, what topology?
#   - output in file/plot in some clever way
#       - maybe as chromosome "bands"?
#       - maybe as % syntenic to each query? => clustering/rooted tree
#       - use plotsr?

def find_pansyn(fins, ref="a", qry="b"):
    """
    TODO docs
    :param:
    :return:
    """
    # columns to look for as start/end positions
    refchr = ref + "chr"
    refstart = ref + "start"
    refend = ref + "end"
    qrychr = qry + "chr"
    qrystart = qry + "start"
    qryend = qry + "end"

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

        while True:
            try: # python iterators suck, so this loop is entirely try-catch'ed
                # this may be very slow, TODO benchmark

                # this uses lexicalic comparisons on strings and is most likely fairly slow
                # TODO possible improvement: store chromosomes as unsigned byte (cython)
                if rsyn[refchr] > lsyn[refchr]:
                    lsyn = next(liter)[1]
                    continue
                if lsyn[refchr] > rsyn[refchr]:
                    rsyn = next(riter)[1]
                    continue
                
                # determine overlap region
                ovstart = max(rsyn[refstart], lsyn[refstart])
                ovend = min(rsyn[refend], lsyn[refend])

                if ovstart < ovend: # there is valid overlap
                    ret = pd.concat([ret, pd.DataFrame(data=[[lsyn[refchr], ovstart, ovend]], columns=left.columns)])
                    lsyn = next(liter)[1]
                    rsyn = next(riter)[1]
                # no valid overlap; one must be in front of the other
                elif lsyn[refstart] > rsyn[refend]: # left is after right
                    rsyn = next(riter)[1]
                else: # right is after left
                    lsyn = next(liter)[1]

            except StopIteration: # nothing more to match
                break

        del riter
        del liter
        return ret


    pansyns = functools.reduce(intersect_syns, syns)

    return pansyns


