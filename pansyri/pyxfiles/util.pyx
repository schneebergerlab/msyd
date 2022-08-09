#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
import multiprocessing
import pandas as pd
from collections import deque
import os

from pansyri.cigar import Cigar
import pansyri.ingest as ingest
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

def extract_regions(fin, ref='a', anns=['SYN'], reforg='ref', qryorg='qry'):
    """
    Given a syri output file, extract all regions matching a given annotation.
    """
    # columns to look for as start/end positions
    refchr = ref + "chr"
    refhaplo = "NaN"
    refstart = ref + "start"
    refend = ref + "end"

    qry = 'b' if ref == 'a' else 'a' # these seem to be the only two values in syri output
    qrychr = qry + "chr"
    qryhaplo = "NaN"
    qrystart = qry + "start"
    qryend = qry + "end"

    buf = deque()
    raw, chr_mapping = ingest.readsyriout(fin) #TODO? handle chr_mapping
    raw = pd.concat([raw.loc[raw['type'] == ann] for ann in anns])
    # if implementing filtering later, filter here

    for row in raw.iterrows():
        row = row[1]
        buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
            Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
            ])

    return pd.DataFrame(data=list(buf), columns=[reforg, qryorg])

def extract_regions_to_list(fins, **kwargs):
    """
    `extract_regions`, but for processing a list of inputs
    """
    return [extract_regions(fin, **kwargs,\
            reforg=fin.split('/')[-1].split('_')[0],\
            qryorg=fin.split('/')[-1].split('_')[-1].split('syri')[0])\
            for fin in fins]

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

