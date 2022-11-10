#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
import multiprocessing
import pandas as pd
from collections import deque
import os
import io
import logging

import pansyri.pansyn as pansyn
from pansyri.classes.cigar import Cigar
from pansyri.classes.coords import Pansyn, Range

logger = logging.getLogger(__name__)

# copied from https://stackoverflow.com/questions/50878960/parallelize-pythons-reduce-command
# doesn't seem to be very fast?
def parallel_reduce(reduceFunc, l, numCPUs):
    if numCPUs == 1 or len(l) <= 3:
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

def parse_input_tsv_path(path):
    """DEPRECATED
    Convenience wrapper to call parse_input_tsv with a path.
    """
    with open(path, 'r') as fin:
        return parse_input_tsv(fin)

def parse_input_tsv(fin):
    """
    Takes a file containing the input alignments/syri files and processes it for find_multisyn.
    Anything after a # is ignored. Lines starting with # are skipped.
    :params: `os.PathLike`, `str` or a TextIO object containing the paths of the input alignment and syri files in tsv format.
    Will be consumed by this function!
    :returns: a tuple of two lists containing the paths of the alignment and syri files.
    """
    if isinstance(fin, (str, os.PathLike)):
        fin = open(fin, 'rt')
    elif not isinstance(fin, io.TextIOBase):
        raise ValueError(f"{fin} is not a path-like or file-like object!")

    syris = deque()     # Lists are too slow appending, using deque instead
    alns = deque()
    for line in fin:
        if line[0] == '#' or line.strip() == '':
            continue

        val = line.strip().split('#')[0].split('\t')
        if len(val) > 2:
            logger.error(f"invalid entry in {fin.name}. Skipping line: {line}")
            continue
        # Check that the files are accessible
        if not os.path.isfile(val[0]):
            raise FileNotFoundError(f"Cannot find file at {val[0]}. Exiting")
        if not os.path.isfile(val[1]):
            raise FileNotFoundError(f"Cannot find file at {val[1]}. Exiting")

        alns.append(val[0].strip())
        syris.append(val[1].strip())

    del fin # maybe just close instead?

    return (syris, alns)


# set of utility funcitons for calling a few preset configurations of find_multisyn using either a list of syri/aln files directly or a tsv containing this information
# For more information, see the find_multisyn docstring
def coresyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), only_core=True, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*parse_input_tsv(path), only_core=False, **kwargs)
def coresyn_from_lists(syns, alns, **kwargs):
    return pansyn.find_multisyn(syns, alns, only_core=True, **kwargs)
def crosssyn_from_lists(syns, alns, **kwargs):
    return pansyn.find_multisyn(syns, alns, only_core=False, **kwargs)

def get_orgs_from_df(df):
    """Small utility function to get all organism from a DataFrame of `Pansyn` objects.
    :param df: A `DataFrame` containing `Pansyn` objects.
    :param_type df: `pandas.DataFrame`
    :returns: A `set` containing all of the organisms in a DataFrame of `pansyn` objects.
    """
    return functools.reduce(lambda x, y: x.union(y), map(lambda x: set(x[1][0].ranges_dict.keys()), df.iterrows()))


#TODO implement function to direcly filter multisyn out for degrees

def get_len(df):
    return sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], df.iterrows())))

def tabularize_lens(df):
    maxdegree = max(map(lambda x: x[1][0].get_degree(), df.iterrows()))
    return [sum(map(lambda x: len(x.ref), filter(lambda x: x.get_degree() == i, map(lambda x: x[1][0], df.iterrows())))) for i in range(maxdegree)]

def get_stats(df):
    """Utility function to output some stats for a df containing computed pansyn objects.
    Calls get_len and tabularize_lens, prettyprints their output.
    """
    tot_len = get_len(df)
    tablens = tabularize_lens(df)
    ret = "Total length: " + str(tot_len) + f" ({len(df)} Regions)" + "\nTotal length by degree:\t" +\
            "\t".join([str(i) for i in range(len(tablens))]) +"\n" +\
            "\t".join([str(i) for i in tablens]) + "\n"
    return ret

def get_call_stats(syns, alns, **kwargs):
    """Utility function to call multisyn in a dataset and immediately compute the statistics using get_stats
    """
    df = pansyn.find_multisyn(syns, alns, **kwargs)
    return get_stats(df)


def eval_combinations(syns, alns, cores=1):
    """Perform get_call_stats for all possible combinations of only_core and SYNAL
    """
    ret = ""
    for only_core, SYNAL in [(False, False), (False, True), (True, False), (True, True)]:
        ret += "core" if only_core else "all"
        ret += " pansynteny, "
        ret += "exact" if SYNAL else "approximate"
        ret += ":\n"
        ret += get_call_stats(syns, alns, only_core=only_core, SYNAL=SYNAL)
    return ret


def filter_multisyn_df(df, rng, only_contained=False):
    """Misc function for filtering a DF produced by find_multisyn for a certain range.
    Only the position on the reference is taken into account.
    Only the chromosome, start and end of the supplied `Range` are used, org and chromosome information is discarded.

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
    #print(inds)
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



