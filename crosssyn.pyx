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

