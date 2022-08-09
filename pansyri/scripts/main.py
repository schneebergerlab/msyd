#!/usr/bin/python3
# in python, probably not worth cythonizing

from pansyri.convenience import *
from pansyri.ordering import *
import pansyri.util as util

import logging
import logging.config

#import argparse as ap
import sys

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

def main(argv):
    if argv[0] == 'order':
        scores = [syn_score, sv_score, len_correct(syn_score), len_correct(sv_score)]
        orgs = argv[1:]
        print("syn_score")
        print(order_greedy(orgs, score_fn=syn_score))
        print("syn_score, len-corrected")
        print(order_greedy(orgs, score_fn=len_correct(syn_score)))
        print("sv_score")
        print(order_greedy(orgs, score_fn=sv_score, maximize=False))
        print("sv_score, len-corrected")
        print(order_greedy(orgs, score_fn=len_correct(sv_score), maximize=False))
        sys.exit()

    syns, alns = util.parse_input_tsv(argv[0])
    cores = int(argv[1]) if len(argv) >= 4 else 1

    if argv[0] == 'len':
        length_compare(syns, alns, cores=cores)
    else:
        eval_combinations(syns, alns, cores=cores)

