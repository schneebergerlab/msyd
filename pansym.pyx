#!/usr/bin/python3
# python right now, convert to cython later

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""
import functools
import pandas as pd
import sys
import ingest


#TODO actually find pansynteny:
#   - use pairwise to ref AND SyRI output (TODO loading)
#   - identify syntenic regions to ref, possibly reusing code from syri(synsearch.pyx, ll. 500-550)
#       - find regions sharing similar locations on the A genome, matching on aStart?
#       - sorted join type algorithm?
#       - lifting pansyntenic regions? if lifting, what topology?
#   - FIRST see if simple intersectionworks, using np.intersect1d
#   - output in file/plot in some clever way
#       - maybe as chromosome "bands"?
#       - maybe as % syntenic to each query? => clustering/rooted tree

#TODO allow loading multiple bam files in main.py

def find_pansyn(fins, ref="a", qry="b"):
    """
    TODO docs
    :param:
    :return:
    """
    # columns to look for as start/end positions
    refstart = ref + "start"
    refend = ref + "end"
    qrystart = qry + "start"
    qryend = qry + "end"

    def extract_syn_regions(df):
        return df.loc[df['type']=='SYN'][["achr", "astart", "aend"]]

    syns = [extract_syn_regions(ingest.readsyriout(fin)[0]) for fin in fins]
    
    pansyns = functools.reduce(pd.merge, syns)

    return len(pansyns)


# helper function to extract only large, syntenic regions suitable for analysis from a DF

print(find_pansyn(sys.argv[1:]))
