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

import argparse
import sys

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

def main(argv):

    #if argv[0] == 'impute':
    #    print(imputation.impute_strings(argv[1], argv[2]))
    #    sys.exit()

    parser = argparse.ArgumentParser(description="Pansyri is a pansynteny and rearrangement identifier.")
    parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Defaults to 4.", type=int, default=4)
    #argparse.FileType('r')
    parser.add_argument("-i", dest='infile', required=True, type=str, help="The .tsv file to read SyRI output and alignment files in from. For more details, see the Readme.")

    inputmode = parser.add_mutually_exclusive_group(required=True)
    inputmode.add_argument("--combinations", "--comb", "--combs", dest='call', action='store_const', const=combinations, help="evaluate all combinations of parameters for pansyn identification.")
    inputmode.add_argument("--lengths", "--len", dest='call', action='store_const', const=lengths, help="Print a table containing the total combined lengths of all pansyntenic regions with a certain degree.")
    inputmode.add_argument("--crossprint", dest='call', action='store_const', const=crossprint, help="Print a DataFrame containing all cross/pansyntenic regions")
    inputmode.add_argument("--coreprint", dest='call', action='store_const', const=coreprint, help="Print a DataFrame containing the core syntenic regions.")
    inputmode.add_argument("--plot", dest='call', action='store_const', const=plot, help="Plotting call, for debugging and validation purposes.")
    inputmode.add_argument("--order", dest='call', action='store_const', const=order, help="Determine the optimal ordering of the supplied genomes for plotting.") #TODO have this be a subparser, accept arguments like which score to use etc

    #TODO separate this into more modular parsing, core vs cross and evaluation separate from that


    args = parser.parse_args(argv)
    args.syns, args.alns = util.parse_input_tsv(args.infile)
    args.call(args)

# call the plotsr ordering functionality on a set of organisms described in the .tsv
def order(args):
    print(ordering.order(args.syns, args.alns))

# compares the output of all four possible values of detect_crosssyn and SYNAL when calling find_multisyn, tabularizes by length
def combinations(args):
    util.eval_combinations(args.syns, args.alns, cores=args.cores)

# do what is done in eval_combinations for every syn/aln file list produced by removing from the end
def lengths(args):
    util.length_compare(args.syns, args.alns, cores=args.cores)

# just print the called cross synteny 
def crossprint(args):
    df = util.crosssyn_from_lists(args.syns, args.alns, cores=args.cores)
    print(df.to_string())

# just print the called core synteny 
def coreprint(args):
    df = util.coresyn_from_lists(args.syns, args.alns, cores=args.cores)
    print(df.to_string())

def plot(args):
    df = util.crosssyn_from_lists(args.syns, args.alns, SYNAL=True)
    cols = ['ref', 'chr'] + list(util.get_orgs_from_df(df))

    def pstolendf(x):
        ret = {k:len(v) for k, v in x.ranges_dict.items()}
        ret['ref'] = len(x.ref)
        ret['chr'] = x.ref.chr # they are all on the same chromosome, so this doesn't matter
        ret = pd.DataFrame(data=ret, columns=cols, index=[0]).fillna(0)
        return ret

    lensdf = pd.concat([pstolendf(x[1][0]) for x in df.iterrows()])
    violating = 0
    violating_length = 0
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
            violating += 1
            violating_length += lens[0]
            logger.error(f"Violates condition: {lens}")
    logger.error(violating)
    sys.exit()
    print(lensdf.to_string(index=False))



