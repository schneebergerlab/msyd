#!/usr/bin/python3
# -*- coding: utf-8 -*-
# dist: language = c++
# cython: language_level = 3

import functools
import multiprocessing
import pandas as pd
from collections import deque
import os

import pansyri.pansyn as pansyn
from pansyri.cigar import Cigar
from pansyri.coords import Pansyn, Range

# copied from https://stackoverflow.com/questions/50878960/parallelize-pythons-reduce-command
# doesn't seem to be very fast?
def parallel_reduce(reduceFunc, l, numCPUs):
    if numCPUs == 1 or len(l) <= 100:
            returnVal = functools.reduce(reduceFunc, l[1:], l[0])
            return returnVal

    parent1, child1 = multiprocessing.Pipe()
    parent2, child2 = multiprocessing.Pipe()
    p1 = multiprocessing.Process(target=parallel_reduce, args=(reduceFunc, l[:len(l) // 2], numCPUs // 2, child1, ) )
    p2 = multiprocessing.Process(target=parallel_reduce, args=(reduceFunc, l[len(l) // 2:], numCPUs // 2 + numCPUs%2, child2, ) )
    p1.start()
    p2.start()
    leftReturn, rightReturn = parent1.recv(), parent2.recv()
    p1.join()
    p2.join()
    returnVal = reduceFunc(leftReturn, rightReturn)
    return returnVal

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
            if line[0] == '#' or line.strip() == '':
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


def coresyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), detect_crosssyn=False, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), detect_crosssyn=True, **kwargs)

maxdegree = 10
len_getter = lambda df: sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], df.iterrows())))
len_tabularizer = lambda df: [sum(map(lambda x: len(x.ref), filter(lambda x: x.get_degree() == i, map(lambda x: x[1][0], df.iterrows())))) for i in range(maxdegree)]

def eval_combinations(syns, alns, cores=1):
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


def length_compare(syns, alns, cores=1):
    syns, alns = list(syns), list(alns)
    for _ in range(len(syns)):
        eval_combinations(syns, alns)
        syns = syns[1:]
        alns = alns[1:]
