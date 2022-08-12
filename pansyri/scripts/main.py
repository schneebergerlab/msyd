#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyri.util as util
import pansyri.ordering as ordering
import pansyri.imputation as imputation
from pansyri.classes.coords import Range

import logging
import logging.config

#import argparse as ap
import sys

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

def main(argv):

    if argv[0] == 'impute':
        print(imputation.impute_strings(argv[1], argv[2]))
        sys.exit()


    syns, alns = util.parse_input_tsv(argv[0])
    cores = int(argv[2]) if len(argv) >= 4 else 1

    # call the plotsr ordering functionality on a set of organisms described in the .tsv
    if argv[1] == 'order': # example call in the ampril dataset folder: pansyri pansr.tsv order 8
        print(ordering.order(syns, alns))
    # compares the output of all four possible values of detect_crosssyn and SYNAL when calling find_multisyn, tabularizes by length
    elif argv[1] == 'comb':# example call in the ampril dataset folder: pansyri pansr.tsv order
        util.eval_combinations(syns, alns, cores=cores)
    # do what is done in eval_combinations for every syn/aln file list produced by removing from the end
    elif argv[1] == 'len':# example call in the ampril dataset folder: pansyri pansr.tsv len 8
        util.length_compare(syns, alns, cores=cores)
    # just print the called cross synteny 
    elif argv[1] == 'print':
        df = util.crosssyn_from_lists(syns, alns, cores=cores)
        print(df)
        print(util.filter_multisyn_df(df, Range(None, 'Chr5', 'NaN', 26000000, 27000000)))

