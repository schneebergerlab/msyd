#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyn
#import Cigar from cigar

import logging
import logging.config

from collections import deque
import argparse as ap
import os
import sys
import numpy as np

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

def coresyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), detect_crosssyn=False, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), detect_crosssyn=True, **kwargs)


def parse_input_tsv(path):
    """
    Takes a file containing the input alignments/syri files and processes it for find_multisyn.
    Anything after a # is ignored. Lines starting with # are skipped.
    :params: path to a file containing the paths of the input alignment and syri files in tsv format
    :returns: a tuple of two lists containing the paths of the alignment and syri files.
    """
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


#import math

def syn_score(syri):
    """
    Defines a distance metric from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    syns = pansyn.extract_regions(syri, ann='SYN')
    return sum(map(lambda x: len(x[1][0]), syns.iterrows()))

def order_plotsr_greedy(orgs, score_fn=syn_score, filename_mapper=lambda x, y: x+'_'+y+"syri.out"):
    """
    A simple, greedy algorithm ordering a list of organisms while trying to maximize the similarity score between each organism and the next one.
    :params:
        orgs is a sequence of organism/filenames
        score sets the similarity score to use, by default `syn_score`. The scoring function needs to accept the filename of a syri.out file as output.
        filename_mapper turns the names of two organisms into the filename of the corresponding syri output.
    :returns: a list containing the elements of orgs ordered according to the algorithm.
    """
    orgs = set(orgs)

    cur = list(orgs)[0] # arbitrarily choose start point
    order = [cur]
    orgs.remove(cur)

    while orgs:
        # find the next organism with maximal similarity score to this one
        max_score = 0# math.inf
        for org in orgs:
            score = score_fn(filename_mapper(cur, org))
            if score > max_score:
                max_score = score
                cur = org

        orgs.remove(cur)
        order.append(cur)

    return order

if sys.argv[1] == 'order':
    orgs = sys.argv[2:]
    print(order_plotsr_greedy(orgs))
    sys.exit()

df1 = coresyn_from_tsv(sys.argv[1], cores=int(sys.argv[2]) if len(sys.argv) >= 3 else 1, sort=True)
#print(df1.to_string())
print("regions:", len(df1))
print("total lengths:", sum(map(lambda x: len(x[1][0].ref), df1.iterrows())))



