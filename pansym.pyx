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
    #   - use pairwise to ref AND SyRI output (TODO loading)
    #   - identify syntenic regions to ref, possibly reusing code from syri(synsearch.pyx, ll. 500-550)
    #       - find regions sharing similar locations on the A genome, matching on aStart?
    #       - sorted join type algorithm?
    #       - lifting pansyntenic regions? if lifting, what topology?
    #   - FIRST see if simple intersectionworks, using np.intersect1d
    #   - output in file/plot in some clever way
    #       - maybe as chromosome "bands"?
    #       - maybe as % syntenic to each query? => clustering/rooted tree
    pass



# pasted from plotsr, parsing syri output
VARS = ['SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
def readsyriout(f):
    from pandas import DataFrame
    import numpy as np
    from collections import deque, OrderedDict
    import logging
    # Reads syri.out. Select: achr, astart, aend, bchr, bstart, bend, srtype
    logger = logging.getLogger("readsyriout")
    syri_regs = deque()
    skipvartype = ['CPG', 'CPL', 'DEL', 'DUPAL', 'HDR', 'INS', 'INVAL', 'INVDPAL', 'INVTRAL', 'NOTAL', 'SNP', 'SYNAL', 'TDM', 'TRANSAL']
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[10] in VARS:
                syri_regs.append(l)
            else:
                if l[10] not in skipvartype:
                    skipvartype.append(l[10])
                    logger.warning("{} is not a valid annotation for alignments in file {}. Alignments should belong to the following classes {}. Skipping alignment.".format(l[10], f, VARS))

    try:
        df = DataFrame(list(syri_regs))[[0, 1, 2, 5, 6, 7, 10]]
    except KeyError:
        raise ImportError("Incomplete input file {}, syri.out file should have 11 columns.".format(f))
    df[[0, 5, 10]] = df[[0, 5, 10]].astype(str)
    try:
        df[[1, 2, 6, 7]] = df[[1, 2, 6, 7]].astype(int)
    except ValueError:
        raise ValueError("Non-numerical values used as genome coordinates in {}. Exiting".format(f))
    # chr ID map
    chrid = []
    chrid_dict = OrderedDict()
    for i in np.unique(df[0]):
        chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]))
        chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]
    df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']
    return df, chrid_dict
# END

