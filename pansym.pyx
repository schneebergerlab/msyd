#!/usr/bin/python3
# python right now, convert to cython later

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""
import functools
import pandas as pd
import numpy as np
import sys
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
    def intersect_syns(left, right, tolerance=10):
        #return pd.merge(left, right) # equivalent to old strict intersection
        
        # they should already be sorted -- required for the rest of the algorithm
        #left.sort_values([refchr, refstart, refend])
        #right.sort_values([refchr, refstart, refend])

        ret = pd.DataFrame(columns = left.columns) # or delete from left? might be faster
        ret = []

        #print(left)
        #print(right)
        riter = right.iterrows()
        for lsyn in left.iterrows():
            try:
                rsyn = next(riter)
            except StopIteration:
                break
            # for some reason, this iterates over tuples with the index
            lsyn = lsyn[1]
            rsyn = rsyn[1]

            if lsyn[refchr] == rsyn[refchr]:
                # left subregion of right
                if lsyn[refstart] > rsyn[refstart] - tolerance and \
                    lsyn[refend] < rsyn[refend] + tolerance:
                        #ret = pd.concat([ret, lsyn])
                        ret.append(lsyn)
                # right subregion of left
                elif rsyn[refstart] > lsyn[refstart] - tolerance and \
                        rsyn[refend] < lsyn[refend] + tolerance:
                        ret.append(rsyn)

        del riter
        return ret


    pansyns = functools.reduce(intersect_syns, syns)

    return len(pansyns)


print(find_pansyn(sys.argv[1:]))
