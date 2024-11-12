#%%cython
#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import sys

import pandas as pd
import numpy as np
import logging
from collections import deque, defaultdict
from multiprocessing import Pool

import msyd.intersection as intersection
from msyd.multisyn import Multisyn, Private
from msyd.coords import Range

cdef MIN_PRIV_THRESH = intersection.MIN_SYN_THRESH


logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)

cpdef complement_dict(syndict):
    """
    Interface to call get_privates on a dict containing multiple multisyn dfs.
    """
    #TODO implement
    pass

cdef complement(multisyns):
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
        rng = msyn.ref
        if rng.start - cov >= MIN_PRIV_THRESH:
            ret.append(Private(Range(rng.org, rng.chr, cov, rng.start - 1)))
        cov = rng.end

    # will ignore private regions at the end of the chromosome;
    # probably makes sense, though

    return pd.DataFrame(ret)


cdef intersect_private(l, r):
    """
    Takes two DataFrames of private regions on the same reference, returns a DF containing regions private in both.
    """
    return find_overlaps(l, r, only_core=True, trim=False, allow_private=True)




#TODO integrate this into main.py/util.py
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

