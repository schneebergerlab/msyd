#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
import copy
import pandas as pd
import numpy as np
import ingest
import util
from cigar import Cigar
import syntools
from syntools import Range, Pansyn
from collections import deque


"""
Notes/TODOs

– properly ctype, including range/positions
then write analogous cross synteny class

use find_multisynteny as entry method, should work as-is

intersect will be broadly similar, but keep & increase count instead of dropping.
Do intersect-add and drop in one step, faster and less messy
=> Benchmark, may be very slow
"""

def combine(l: Pansyn, r: Pansyn, keep_overlap=False):
    #TODO implement
    # use mostly the same overlap calculation
    # but replacing old objects should not be necessary? => discuss

    # compute overlap
    ovstart = max(l.ref.start, r.ref.start)
    ovend = min(l.ref.end, r.ref.end)
    assert(ovstart < ovend, f"ERROR: no overlap found between {l} and {r}")

    leftest = l if l.ref.start < r.ref.start else r
    rightest = l if l.ref.end > r.ref.end else r

    ret = set() # use a set to automatically remove duplicates

    if keep_overlap:
        ret.add(leftest)
    else:
        ret.add(leftest.drop(0, leftest.ref.end - ovstart))
    
    # core synteny
    ret.add(l.drop(ovstart - l.ref.start, l.ref.end - ovend) + r.drop(ovstart - r.ref.start, r.ref.end - ovend))

    if keep_overlap:
        ret.add(rightest)
    else:
        ret.add(rightest.drop(0, rightest.ref.end - ovstart))

    return sorted(ret)



# 



def intersect_crosssyns(left, right):
    pass
