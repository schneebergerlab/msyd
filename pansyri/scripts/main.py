#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyri.util as util
import pansyri.ordering as ordering
import pansyri.imputation as imputation

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

    if argv[1] == 'order':
        print(ordering.order(syns, alns))
    elif argv[1] == 'len':
        util.length_compare(syns, alns, cores=cores)
    elif argv[1] == 'comb':
        util.eval_combinations(syns, alns, cores=cores)
    elif argv[1]:
        print(util.crosssyn_from_lists(syns, alns, cores=cores))
