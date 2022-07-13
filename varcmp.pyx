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
    
    :params:
        `syri`: dataframe containing the information from a syri file, as parsed by `ingest.readsyriout`
        `rngs`: a List of `Range` objects to extract snps from.
    :returns:
        A dictionary mapping a `Pansyn` to a dictionary containing the corresponding SNPs identified on it for each organism that `Pansyn` has a position in.
    """
    pass

# idea from manish: check equality between SVs by neighbourhood pansyntenic regions
def get_pansyn_neighbours(rng: Range, pansyns, ref=True):
    """
    Given a `Range` and a list of pansyntenic regions (as `Pansyn`), returns the regions neighbouring or overlapping with that region on the organism specified in the `Range`.
    :returns: a list of `Pansyn` objects
    """
    pass

def cmp_neighbourhoods(l, r):
    """
    Given two neighbourhoods calculated by `get_pansyn_neighbours`, tries to determine if they are the same.
    :returns: True/False, optionally maybe output distance value
    """
    pass
