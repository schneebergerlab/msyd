#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyri.util as util
import pansyri.ordering as ordering
import pansyri.ingest as ingest
import pansyri.imputation as imputation
import pansyri.pansyn as pansyn
from pansyri.classes.coords import Range

import logging
import logging.config

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

import pandas as pd

import argparse
import sys

"""
This file serves as the main entrypoint for the pansyri CLI.
"""

def main(argv):

    #if argv[0] == 'impute':
    #    print(imputation.impute_strings(argv[1], argv[2]))
    #    sys.exit()

    ## Argument parsing redo
    # Idea:
    # – different modes, selectablel by first positional argument
    # – call, filter, statistics, order, lengths mode?
    # – have subparsers for each, mb share common arguments like -c and -i/-o


    parser = argparse.ArgumentParser(description="Pansyri is a pansynteny and rearrangement identifier.")
    # core arguments
    parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Defaults to 4.", type=int, default=4)
    parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="The .tsv file to read SyRI output and alignment files in from. For more details, see the Readme.")
    parser.add_argument("-o", dest='outfile', default='-', type=argparse.FileType('wt'), help="Where to store disk output. Defaults to stdout (specified with \"-\").")

    subparsers = parser.add_subparsers(description="See also pansyri [subparser] -h:") # title/description?
    call_parser = subparsers.add_parser("call", description="Call Pansynteny and write to disk in VCF or PFF format.")
    order_parser = subparsers.add_parser("order", description="Determine the optimal ordering of the supplied genomes for plotting.")
    filter_parser = subparsers.add_parser("filter", description="")
    combinations_parser = subparsers.add_parser("combinations", description="Evaluate all combinations of --syn and --core in the same manner as the --stats parameter.")
    statistics_parser = subparsers.add_parser("statistics", description="Print some statistics of the called pansynteny.")
    plot_parser = subparsers.add_parser("plot", description="Prints a lengths df for plotting to stdout. Can be piped to a file and plotted with tests/plot.R .")


    call_parser.add_argument("--core", dest='core', action='store_const', const=True, default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    call_parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL regions, yields more contiguous regions and faster runtime, but calls may not be exact to the base level.")
    call_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified")
    call_parser.add_argument("-v", "--vcf", dest='vcf', action='store_const', const=False, default=True, help="Store as .vcf file. Discards CIGAR strings!")

    statistics_parser.add_argument("--vcf", dest='invcf', required=False, type=argparse.FileType('r'), help="The .vcf file to filter and write to -o. Only necessary when calling --filter-vcf.")

    args = parser.parse_args(argv)
    args.syns, args.alns = util.parse_input_tsv(args.infile)
    args.df = pansyn.find_multisyn(args.syns, args.alns, only_core=args.core, SYNAL=args.SYNAL)
    if args.call:
        args.call(args)
    else:
        util.to_pff(args.df.head(args.limit) if args.limit != 0 else args.df, args.outfile, save_cigars=args.cigars)



# call the plotsr ordering functionality on a set of organisms described in the .tsv
def order(args): # TODO have subparsing for other ordering/scores/etc
    print(ordering.order_hierarchical(args.df, orgs=None, score_fn=ordering.syn_score))

# wrapper around the util fucntion
def stats(args):
    print(util.get_stats(args.df))

# compares the output of all four possible values of detect_crosssyn and SYNAL when calling find_multisyn, tabularizes by length
def combinations(args):
    print(util.eval_combinations(args.syns, args.alns, cores=args.cores))

# do what is done in eval_combinations for every syn/aln file list produced by removing from the end
def lengths(args):
    util.length_compare(args.syns, args.alns, cores=args.cores)

def filter_vcf(args):
    if not args.invcf:
        logger.error("No vcf to filter specified!")
        raise ValueError("No vcf to filter specified!")
    # close the incoming handles to avoid double-opening the files
    args.invcf.close()
    args.outfile.close()
    ingest.extract_syntenic_from_vcf(args.df, args.invcf, args.outfile.name)


def plot(args):
    df = args.df
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
    print(lensdf.to_string(index=False))



