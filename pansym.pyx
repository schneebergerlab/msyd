#!/usr/bin/python3
# python right now, convert to cython later

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""

def find_pansyn(coords):
    #TODO actually find pansynteny:
    #   - use pairwise to ref or multiple alignments?
    #       => start with pairwise to ref
    #   - identify syntenic regions to ref, possibly reusing code from syri(synsearch.pyx, ll. 500-550)
    #   - output in file/plot in some clever way
    #       - maybe as chromosome "bands"?
    #       - maybe as % syntenic to each query? => clustering/rooted tree
