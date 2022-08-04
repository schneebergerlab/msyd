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


maxdegree = 10
len_getter = lambda df: sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], df.iterrows())))
len_tabularizer = lambda df: [sum(map(lambda x: len(x.ref), filter(lambda x: x.get_degree() == i, map(lambda x: x[1][0], df.iterrows())))) for i in range(maxdegree)]


def length_compare():
    syns, alns = parse_input_tsv(sys.argv[2])
    syns, alns = list(syns), list(alns)
    cores = int(sys.argv[3]) if len(sys.argv) >= 4 else 1
    for _ in range(len(syns)):
        cores_syn = pansyn.find_multisyn(syns, alns, detect_crosssyn=False, sort=True, SYNAL=False, cores=cores)
        cores_synal = pansyn.find_multisyn(syns, alns, detect_crosssyn=False, sort=True, SYNAL=True, cores=cores)

        cross_syn = pansyn.find_multisyn(syns, alns, detect_crosssyn=True, sort=True, SYNAL=False, cores=cores)
        cross_synal = pansyn.find_multisyn(syns, alns, detect_crosssyn=True, sort=True, SYNAL=True, cores=cores)

        print("\nComparing", syns)
        
        print(f"cores_syn:\tRegions: {len(cores_syn)}\tTotal length:{len_getter(cores_syn)}")
        print(f"cores_synal:\tRegions: {len(cores_synal)}\tTotal length:{len_getter(cores_synal)}")
        print("")
        print(f"cross_syn:\tRegions: {len(cross_syn)}\tTotal length:{len_getter(cross_syn)}")
        print("degree:", '\t\t'.join([str(i) for i in range(maxdegree)]), sep='\t')
        print("count:", '\t'.join([str(i) for i in len_tabularizer(cross_syn)]), sep='\t')
        print(f"cross_synal:\tRegions: {len(cross_synal)}\tTotal length:{len_getter(cross_synal)}")
        print("degree:", '\t\t'.join([str(i) for i in range(maxdegree)]), sep='\t')
        print("count:", '\t'.join([str(i) for i in len_tabularizer(cross_synal)]), sep='\t')

        syns = syns[1:]
        alns = alns[1:]

    sys.exit()

if sys.argv[1] == 'order':
    orgs = sys.argv[2:]
    print(order_plotsr_greedy(orgs))
    sys.exit()

if sys.argv[1] == 'len':
    length_compare()

df1 = coresyn_from_tsv(sys.argv[1], cores=int(sys.argv[2]) if len(sys.argv) >= 3 else 1, sort=True, SYNAL=False)
df2 = crosssyn_from_tsv(sys.argv[1], cores=int(sys.argv[2]) if len(sys.argv) >= 3 else 1, sort=True, SYNAL=False)
#print(df1.to_string())
print(df1)
print(df2)
print("regions:", len(df1))
print("total coresyn length:", sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], df1.iterrows()))))



