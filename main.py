#!/usr/bin/python3
# in python, probably not worth cythonizing

import ingest
import coresyn

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



df1 = syntools.coresyn_from_tsv(sys.argv[1], cores=int(sys.argv[2]) if len(sys.argv) >= 3 else 1)
#print(df1.to_string())
print("regions:", len(df1))
print("total lengths:", sum(map(lambda x: x[1][0].end-x[1][0].start, df1.iterrows())))
