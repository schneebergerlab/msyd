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
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

def main(argv):

    #if argv[0] == 'impute':
    #    print(imputation.impute_strings(argv[1], argv[2]))
    #    sys.exit()

    parser = argparse.ArgumentParser(description="Pansyri is a pansynteny and rearrangement identifier.")
    parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Defaults to 4.", type=int, default=4)
    parser.add_argument("--limit", dest="limit", help="Limits the number of lines that are output to an amount reasonable for printing. Defaults to 0 (no limit)!", type=int, default=0)
    parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="The .tsv file to read SyRI output and alignment files in from. For more details, see the Readme.")
    parser.add_argument("--in-vcf", dest='invcf', required=False, type=argparse.FileType('r'), help="The .vcf file to filter and write to -o. Only necessary when calling --filter-vcf.")
    parser.add_argument("-o", dest='outfile', default='-', type=argparse.FileType('wt'), help="Where to send output to in .pff format. Defaults to stdout (specified with \"-\").\nSending output can be disabled with --discard. This parameter is ignored by --stats, --combinations, --lengths and --order, which will always print to stdout. Will be in PFF format, unless --filter-vcf is specified, in which case output will be in vcf format.")
    parser.add_argument("--core", dest='core', action='store_const', const=True, default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL regions, yields more contiguous regions and faster runtime, but calls may not be exact to the base level.")
    parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified")


    inputmode = parser.add_mutually_exclusive_group(required=False)
    inputmode.add_argument("--statistics", "--stats", dest='call', action='store_const', const=stats, help="Print some statistics about the called multisynteny.")
    inputmode.add_argument("--combinations", "--comb", "--combs", dest='call', action='store_const', const=combinations, help="Evaluate all combinations of --syn and --core in the same manner as the --stats parameter.")
    inputmode.add_argument("--lengths", "--len", "--lens", dest='call', action='store_const', const=lengths, help="Debug call to perform the combinations evaluation (--comb) on subsets of the genomes supplied, generated by iteratively removing query genomes from the top down.")
    inputmode.add_argument("--plot", dest='call', action='store_const', const=plot, help="Prints a lengths df for plotting to stdout. Can be piped to a file and plotted with tests/plot.R .")
    inputmode.add_argument("--order", dest='call', action='store_const', const=order, help="Determine the optimal ordering of the supplied genomes for plotting.") #TODO have this be a subparser, accept arguments like which score to use etc
    inputmode.add_argument("--discard", dest='call', action='store_const', const=lambda x: None, help="Discard the computed output. For benchmarking.")
    inputmode.add_argument("--filter-vcf", dest='call', action='store_const', const=filter_vcf, help="Discard the computed output. For benchmarking.")


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



