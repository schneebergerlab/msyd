#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import numpy as np
from pansyn import Range, Pansyn
import pandas as pd

def get_snps(syri, pansyns):
    """
    
    :params:
        `syri`: dataframe containing the information from a syri file, as parsed by `ingest.readsyriout`
        `pansyns`: a List of `Pansyn` objects to extract snps from.
    :returns:
        A dictionary mapping a `Pansyn` to a dictionary containing the corresponding SNPs identified on it for each organism that `Pansyn` has a position in.
    """
    # idea from manish: check equality between SVs by neighbourhood pansyntenic regions
    pass

