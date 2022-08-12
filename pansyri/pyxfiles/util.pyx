#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
import multiprocessing
import pandas as pd
from collections import deque
import os

import pansyri.pansyn as pansyn
from pansyri.classes.cigar import Cigar
from pansyri.classes.coords import Pansyn, Range

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


# set of utility funcitons for calling a few preset configurations of find_multisyn using either a list of syri/aln files directly or a tsv containing this information
# For more information, see the find_multisyn docstring
def coresyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), detect_crosssyn=False, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), detect_crosssyn=True, **kwargs)
def coresyn_from_lists(syns, alns, **kwargs):
    return pansyn.find_multisyn(syns, alns, detect_crosssyn=False, **kwargs)
def crosssyn_from_lists(syns, alns, **kwargs):
    return pansyn.find_multisyn(syns, alns, detect_crosssyn=True, **kwargs)

def get_orgs_from_df(df):
    """Small utility function to get all organism from a DataFrame of `Pansyn` objects.
    :param df: A `DataFrame` containing `Pansyn` objects.
    :param_type df: `pandas.DataFrame`
    :returns: A `set` containing all of the organisms in a DataFrame of `pansyn` objects.
    """
    return functools.reduce(lambda x, y: x.union(y), map(lambda x: set(x[1][0].ranges_dict.keys()), df.iterrows()))


#TODO implement function to direcly filter multisyn out for degrees

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

def filter_multisyn_df(df, rng, only_contained=False):
    """Misc function for filtering a DF produced by find_multisyn for a certain range.
    Only the position on the reference is taken into account.
    Only the chromosome, start and end of the supplied Range are used, org and chromosome information is discarded.

    :param df: `find_multisyn` `DataFrame` of `Pansyn` objects.
    :type df: `DataFrame[Pansyn]`
    :param rng: `Range` for selecting the `Pansyn` objects.
    :param only_contained: switches between selecting any region intersecting or contained in the specified `Range`.
    :type only_contained: bool
    """
    def filter_fn(x):
        ref = x.ref
        # check if on same chr
        if not rng.chr == ref.chr:
            return False
        # check if contained:
        if rng.start < ref.start < rng.end and rng.start < ref.end < rng.end:
            return True
        # quit if only looking for contained regions
        if only_contained:
            return False
        # check if start or end within rng
        return rng.start < ref.start < rng.end or rng.start < ref.end < rng.end

    inds = df[0].apply(filter_fn)
    print(inds)
    return df.loc[inds]

def filter_multisyn_df_chr(df, chr):
    """Misc function for filtering a DF produced by find_multisyn for a certain chromosome.
    Does essentially the same thing as `filter_multsyn_df`, but only uses chromosome information

    :param df: `find_multisyn` `DataFrame` of `Pansyn` objects.
    :type df: `DataFrame[Pansyn]`
    :param chr: Chromosome to select.
    :type chr: `str`
    """
    return df.loc[df[0].apply(lambda x: chr == x.ref.chr)]

def length_compare(syns, alns, cores=1):
    syns, alns = list(syns), list(alns)
    for _ in range(len(syns)):
        eval_combinations(syns, alns)
        syns = syns[1:]
        alns = alns[1:]



