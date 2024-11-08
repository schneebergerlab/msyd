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


logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)

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

