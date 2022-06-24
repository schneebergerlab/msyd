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

def parse_input_tsv(path):
    """
    Takes a file containing the input alignments/syri files and processes it for coresyn.pyx.
    Anything after a # is ignored. Lines starting with # are skipped.
    :params: path to a file containing the paths of the input alignment and syri files in tsv format
    :returns: a tuple of two lists containing the paths of the alignment and syri files.
    """
    from collections import deque
    import os
    syris = deque()     # Lists are too slow appending, using deque instead
    alns = deque()
    with open(path, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue

            val = line.strip().split('#')[0].split('\t')
            if len(val) > 2:
                print(f"ERROR: invalid entry in {path}. Skipping line: {line}")
                continue
            # Check that the files are accessible
            if not os.path.isfile(val[0]):
                raise FileNotFoundError(f"Cannot find file at {val[0]}. Exiting")
            if not os.path.isfile(val[1]):
                raise FileNotFoundError(f"Cannot find file at {val[1]}. Exiting")

            alns.append(val[0].strip())
            syris.append(val[1].strip())

    return (syris, alns)


remcigar = lambda x: x# x[0] if type(x)==list or type(x)==tuple else x

df1 = coresyn.find_coresyn(*parse_input_tsv(sys.argv[1]), sort=False).apply(lambda x: x.apply(remcigar), cores=int(sys.argv[2]) if len(sys.argv) >= 3 else 1)
#print(df1.to_string())
print("regions:", len(df1))
print("total lengths:", sum(map(lambda x: x[1][0].end-x[1][0].start, df1.iterrows())))
