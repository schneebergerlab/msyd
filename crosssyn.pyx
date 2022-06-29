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
from syn import Range


"""
Notes/TODOs

– properly ctype, including range/positions
then write analogous cross synteny class

use find_multisynteny as entry method, should work as-is

intersect will be broadly similar, but keep & increase count instead of dropping.
Do intersect-add and drop in one step, faster and less messy
=> Benchmark, may be very slow
"""

class Crosssyn:
    """
    TODO docstring
    """
    # ranges, cigars have type List[Range]/List[Cigar], respectively, but cython cannot deal with generic type hints
    def __init__(self, ref:Range, ranges, cigars):
        self.ref = ref # optional if using a reference-free algorithm. NONE CURRENTLY IMPLEMENTED!
        self.ranges = ranges
        self.cigars = cigars # length equal to ranges; optional if using approximate matching
        self.degree = len(self.ranges)

    def __repr__(self):
        return f"Crosssyn({self.ref}, {self.ranges})"

    def __eq__(l, r):
        return l.ref == r.ref and l.ranges == r.ranges and l.cigars == r.cigars
        
    # for now, only sorts on the reference (falling back to the Range comparison operator)
    def __lt__(l, r):
        if not l.ref or not r.ref:
            raise ValueError(f"ERROR comparing {l} with {r}: both need to have a reference!")
        return l.ref < r.ref

    def add(self, rng:Range, cg: Cigar):
        self.ranges.append(rng)
        self.cigars.append(cg)
        self.degree += 1

    def combine(l: Crosssyn, r: Crosssyn, keep_overlap=False):
        #TODO implement
        # return list or iterator instead of new object
        # use mostly the same overlap calculation
        # but replacing old objects should not be necessary? => discuss
        # degree calculation should be handled in constructor

        # compute overlap
        ovstart = max(l.ref.start, r.ref.start)
        ovend = min(l.ref.end, r.ref.end)
        assert(ovstart < ovend, f"ERROR: no overlap found between {l} and {r}")

        rdropstart = ovstart - r.ref.start # 0 if r.ref is maximal, else positive
        ldropstart = ovstart - l.ref.start 

        rdropend = r.ref.end - ovend # 0 if r.ref is minimal, else positive
        ldropend = l.ref.end - ovend

        start, overlap, end = None

        if l.cigars and r.cigars:

        else:
            return sorted([
                    Crosssyn(l.ref.drop(0, l.ref.end

                    ])




def intersect_crosssyns(left, right):
