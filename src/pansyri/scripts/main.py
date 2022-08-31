#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyri.util as util
import pansyri.ordering as ordering
import pansyri.imputation as imputation
from pansyri.classes.coords import Range

import logging
import logging.config

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

import pandas as pd

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
        print(ordering.order(syns, alns, argv[2] if len(argv) >= 3 else None))
    # compares the output of all four possible values of detect_crosssyn and SYNAL when calling find_multisyn, tabularizes by length
    elif argv[1] == 'comb':# example call in the ampril dataset folder: pansyri pansr.tsv order
        util.eval_combinations(syns, alns, cores=cores)
    # do what is done in eval_combinations for every syn/aln file list produced by removing from the end
    elif argv[1] == 'len':# example call in the ampril dataset folder: pansyri pansr.tsv len 8
        util.length_compare(syns, alns, cores=cores)
    # just print the called cross synteny 
    elif argv[1] == 'crossprint':
        df = util.crosssyn_from_lists(syns, alns, cores=cores)
        print(df.to_string())
    elif argv[1] == 'coreprint':
        df = util.coresyn_from_lists(syns, alns, cores=cores)
        print(df.to_string())
    elif argv[1] == 'plot':
        df = util.crosssyn_from_lists(syns, alns, SYNAL=True)
        cols = ['ref', 'chr'] + list(util.get_orgs_from_df(df))

        def pstolendf(x):
            ret = {k:len(v) for k, v in x.ranges_dict.items()}
            ret['ref'] = len(x.ref)
            ret['chr'] = x.ref.chr # they are all on the same chromosome, so this doesn't matter
            ret = pd.DataFrame(data=ret, columns=cols, index=[0]).fillna(0)
            return ret

        lensdf = pd.concat([pstolendf(x[1][0]) for x in df.iterrows()])
        for row in lensdf.loc[:, lensdf.columns != 'chr'].iterrows():
            violates = False
            lens = row[1]
            minlen = lens[0]*0.9
            maxlen = lens[0]*1.1
            for l in lens:
                if l > 0:
                    if l < minlen or l > maxlen:
                        violates = True

            if violates:
                logger.error(f"Violates condition: {lens}")
        sys.exit()
        print(lensdf.to_string(index=False))



