#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import numpy as np
from pansyn import Range, Pansyn
import pandas as pd

def get_snps(syri, rngs):
    """
    A function that finds all SNPs annotated within one of the input region, as a dictionary keyed by organism.
    :params:
        `syri`: dataframe containing the information from a syri file, as parsed by `ingest.readsyriout`
        `rngs`: a List of `Range` objects to extract snps from.
    :returns:
        A dictionary mapping a `Pansyn` to a dictionary containing the corresponding SNPs identified on it for each organism that `Pansyn` has a position in.
    """
    pass

#TODO maybe refactor neighbourhoods out into a class to allow separate saving of overlap/ proper neighbours, comparison operators etc.
# maybe optionally save with variant annotation?

# idea from manish: check equality between SVs by neighbourhood pansyntenic regions
def get_pansyn_neighbours(rng: Range, pansyns, ref=True, overlapping=True):
    """
    Given a `Range` and a list of pansyntenic regions (as `Pansyn`), returns the regions neighbouring or overlapping with that region on the organism specified in the `Range`.
    :returns: a list of `Pansyn` objects
    """
    # neighbouring := nearest in each direction
    if len(pansyns) < 2:
        print("WARN: get_pansyn_neighbours called with just one pansyn, returning as neighbours!")
        return pansyns

    org = rng.org
    ret = []
    prev = None
    paniter = iter(pansyns.ranges)
    prev = next(paniter)
    cur = next(paniter)
    try:
        while cur.ref.end < rng.start if ref else cur.ranges_dict[org].end < rng.start: # maybe TODO refactor to use more elegant polymorphism than ternaries?
            prev = cur
            cur = next(paniter)

        ret.append(prev) # left neighbour
        del prev # no longer need to save previous values
        
        # get overlapping 'neighbours'
        while cur.ref.start < rng.end if ref else cur.ranges_dict[org].start < rng.end:
            if overlapping:
                ret.append(cur)
            cur = next(paniter)

        ret.append(cur) # right neighbour
                
    except StopIteration:
        print("WARN: get_pansyn_neighbours found rng to be right of all pansyntenic regions! Double-Check if calling with the right set of pansyns")

    return ret



def cmp_neighbourhoods(l, r):
    """
    Given two neighbourhoods calculated by `get_pansyn_neighbours`, tries to determine if they are the same.
    Use a function to be able to use more sophisticated comparisons.
    :returns: True/False, optionally maybe output distance value
    """
    pass
