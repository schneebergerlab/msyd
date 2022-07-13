#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyn
import Cigar from cigar

import logging
import logging.config

import argparse as ap
import os
import sys
import numpy as np

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

def coresyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*pansyn.parse_input_tsv(path), detect_crosssyn=False, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*pansyn.parse_input_tsv(path), detect_crosssyn=True, **kwargs)


df1 = coresyn_from_tsv(sys.argv[1], cores=int(sys.argv[2]) if len(sys.argv) >= 3 else 1)
#print(df1.to_string())
print("regions:", len(df1))
print("total lengths:", sum(map(lambda x: x[1][0].ref.end-x[1][0].ref.start, df1.iterrows())))


