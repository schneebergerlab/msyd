#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
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

    ret = []

    # problems: order according to start/end
    # detect when no region should be added

    # preliminary, assuming l ist leftest

    ldrop = len(l.ref) - ldropstart
    ret.append(Pansyn(l.ref.drop(0, ldrop), [rng.drop(0, ldrop) for rng in l.ranges], None))


    return sorted(ret)



# 



def intersect_crosssyns(left, right):
    pass
