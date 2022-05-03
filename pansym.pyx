#!/usr/bin/python3
# python right now, convert to cython later

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""

# coords is a list of coordinate pandas frames as output by the readCoords function.
# A is the reference and B the query genome
#TODO allow loading multiple bam files in main.py

def find_pansyn(coords):
    #TODO actually find pansynteny:
    #   - use pairwise to ref or multiple alignments?
    #       => start with pairwise to ref
    #   - identify syntenic regions to ref, possibly reusing code from syri(synsearch.pyx, ll. 500-550)
    #       - find regions sharing similar locations on the A genome, matching on aStart?
    #       - sorted join type algorithm?
    #       - lifting pansyntenic regions? if lifting, what topology?
    #   - output in file/plot in some clever way
    #       - maybe as chromosome "bands"?
    #       - maybe as % syntenic to each query? => clustering/rooted tree
    pass

