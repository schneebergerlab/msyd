#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pansyri.pansyn as pansyn
import pansyri.util as util
import pansyri.cigar as cigar
from pansyri.cigar import Cigar

def coresyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*util.parse_input_tsv(path), detect_crosssyn=False, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pansyn.find_multisyn(*util.parse_input_tsv(path), detect_crosssyn=True, **kwargs)

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
