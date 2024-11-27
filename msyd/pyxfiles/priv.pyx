#%%cython
#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3


import sys
import logging

from collections import deque, defaultdict
from multiprocessing import Pool

import pandas as pd
import numpy as np

import msyd.intersection as intersection
import msyd.util as util
from msyd.multisyn import Multisyn, Private
from msyd.coords import Range

logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)


cdef int MIN_PRIV_THRESH = intersection.get_min_syn_thresh()



#THOUGHT implement a higher-order function for handling these kinds of things in utils?
cpdef complement_dict(syndict, add=False, cores=1):
    """
    Interface to call complement on a dict containing multiple multisyn DFs.
    :returns: A dict mapping chromosomes to DFs containing the private and/or Multisyn objects on them.
    """
    if cores > 1:
        with Pool(cores) as pool:
            return dict(pool.map(_workaround, [(it[0], it[1], add) for it in syndict.items()]))
    else:
        return dict(map(_workaround, [(it[0], it[1], add) for it in syndict.items()]))

cpdef _workaround(tup): # tup: [chrom, msyn, add]
    # Annoying workaround, because multiprocessing doesn't like lambdas
    return tup[0], complement(tup[1], add=tup[2])

cdef complement(multisyns, add=False):
    """
    Takes a DF of multisyn objects, finds regions on the reference not covered by any of them.
    Assumes the multisyns are sorted by position on the reference.
    This can be thought of as a way to compute the complement of the syntenic regions.

    :returns: a DF containing Private objects, sorted by position on reference
    """
    #THOUGHT: should private inherit from multisyn?
    # could also just use Ranges for it, code in intersection should be easy to copy/adopt

    cdef:
        int cov = 0
        ret = deque()
    
    for _, msyn in multisyns.iterrows():
        msyn = msyn[0]

        # compute gap between previous msyn and this, annotate if large enough
        rng = msyn.ref
        if rng.start - cov >= MIN_PRIV_THRESH:
            # ranges are inclusive on both end and start
            ret.append(Private(Range(rng.org, rng.chr, cov + 1, rng.start - 1)))
        cov = rng.end
        
        # if specified, also add the msyn while maintaining sorting
        if add:
            ret.append(msyn)

    # will ignore private regions at the end of the chromosome;
    # probably hard to explicitly annotate them due to unknown telomere length anyway

    return pd.DataFrame(data=list(ret))

cdef intersect_private(l, r):
    """
    Takes two DataFrames of private regions on the same reference, returns a DF containing regions private in both.
    """
    return intersection.find_overlaps(l, r, only_core=True, trim=False, allow_private=True)

#TODO debug/validate
#TODO make complement work with realignment datastructure
# => will require accounting for mappingtree/other multisyns 
# => do this on the synthetic sequences, before mapping back to real sequence space?
# => spacers might be a problem


# Idea for additional fn
# finds private regions by scanning through the genome for regions not covered by any merasyn
# tracks current position along the genome
# challenge: non-coresyn regions
# => approach: sort merasyns by org, then subtract
# alternatively, use intervaltrees, subtract each merasyn
# maybe move to own file eventually?
# implement after refactoring of data structures
#cpdef find_private(syns, only_private=False):
#    # Finds 
#    pass

