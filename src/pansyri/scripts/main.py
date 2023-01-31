#!/usr/bin/python3
# in python, probably not worth cythonizing

import pansyri.util as util
import pansyri.ordering as ordering
import pansyri.io as io
import pansyri.imputation as imputation
import pansyri.pansyn as pansyn
from pansyri.classes.coords import Range

logger = util.CustomFormatter.getlogger(__name__)

import pandas as pd

import argparse
import sys

"""
This file serves as the main entrypoint for the pansyri CLI.
"""

def main(argv):

    parser = argparse.ArgumentParser(description="Pansyri is a pansynteny and rearrangement identifier.")
    parser.set_defaults(func=None, cores=1)

    subparsers = parser.add_subparsers(description="See also pansyri [subparser] -h:") # title/description?
    # ordering parser
    order_parser = subparsers.add_parser("order", description="Determine the optimal ordering of the supplied genomes for plotting.")
    order_parser.set_defaults(func=order)
    order_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")

    #combinations_parser = subparsers.add_parser("combinations", description="Evaluate all combinations of --syn and --core in the same manner as the --stats parameter.")
    #combinations_parser.set_defaults(func=combinations)
    #combinations_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")

    ## plotting subparser
    #plot_parser = subparsers.add_parser("plot", description="Prints a lengths df for plotting to stdout. Can be piped to a file and plotted with tests/plot.R .")
    #plot_parser.set_defaults(func=plot)
    #plot_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF or VCF file to read pansynteny information from.")

    # Pansyn calling argparser
    call_parser = subparsers.add_parser("call", help="test", description="Call Pansynteny and write to disk in VCF or PFF format.")
    call_parser.set_defaults(func=call)
    call_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="The .tsv file to read SyRI output and alignment files in from. For more details, see the Readme.")
    call_parser.add_argument("-o", dest='pff', type=argparse.FileType('wt'), help="Where to save the output in PFF format (see format.md). At least one of -o and -v must be specified!")
    call_parser.add_argument("-v", "--vcf", dest='vcf', type=argparse.FileType('wt'), help="Where to save the output in VCF format. At least one of -o and -v must be specified!")
    call_parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Defaults to 4.", type=int, default=4)
    call_parser.add_argument("--core", dest='core', action='store_const', const=True, default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    call_parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL regions, yields more contiguous regions and faster runtime, but calls may not be exact to the base level.")
    call_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified")
    call_parser.add_argument("-p", "--print", dest='print', action='store_true', default=False, help="print a subset of the output to stdout, for debugging.")

    # Conversion parser
    convert_parser = subparsers.add_parser("convert", description="Convert between different pansynteny information file formats")
    convert_parser.set_defaults(func=convert)
    convert_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    convert_parser.add_argument("-t", dest='filetype', required=False, type=str, help="File format to convert to.")
    convert_parser.add_argument("-o", dest='outfile', default='-', type=argparse.FileType('wt'), help="Where to store disk output. Defaults to stdout (specified with \"-\").")



    # VCF filtering subparser
    view_parser = subparsers.add_parser("view", help="", description="")
    view_parser.set_defaults(func=view)
    view_parser.add_argument("--vcf", dest='invcf', required=True, type=argparse.FileType('r'), help="The .vcf file to filter and write to -o.")
    view_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    view_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the filtered VCF file. Defaults to stdout (specified with \"-\").")

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)
    else:
        #logger.info("No subcommand specified, priting help message.")
        parser.print_help()

def call(args):
    qrynames, syns, alns = util.parse_input_tsv(args.infile)
    df = pansyn.find_multisyn(qrynames, syns, alns, only_core=args.core, SYNAL=args.SYNAL)
    if args.print:
        logger.info("Printing sample head to STDOUT")
        print(df.head())

    print(util.get_stats(df))

    # save output
    logger.info("Saving pansyn calls to file(s)")
    if not ( args.vcf or args.pff):
        logger.error("At least one output file should be specified! Quitting!")
        sys.exit(-1)

    if args.pff:
        io.save_to_pff(df, args.pff, save_cigars=args.cigars)
    if args.vcf:
        io.save_to_vcf(df, args.vcf, cores=args.cores)

# call the plotsr ordering functionality on a set of organisms described in the .tsv
def order(args):
    df = io.read_pff(args.infile)
    print(ordering.order_hierarchical(df, orgs=None, score_fn=ordering.syn_score))

# compares the output of all four possible values of detect_crosssyn and SYNAL when calling find_multisyn, tabularizes by length
def combinations(args):
    syns, alns = util.parse_input_tsv(args.infile)
    print(util.eval_combinations(syns, alns, cores=args.cores))

# do what is done in eval_combinations for every syn/aln file list produced by removing from the end
def lengths(args):
    syns, alns = util.parse_input_tsv(args.infile)
    util.length_compare(syns, alns, cores=args.cores)

def convert(args):
    df = io.read_pff(args.infile)
    # mb refactor into common output function
    if not args.filetype:
        args.filetype = args.outfile.name.split(".")[-1]
    if args.filetype not in {'vcf', 'pff'}: # add bcf
        logger.error("invalid output file format selected!")
    if args.filetype == 'pff':
        io.save_to_pff(df, args.outfile)
    else:
        io.save_to_vcf(df, args.outfile, cores=args.cores)

def view(args):
    df = io.read_pff(args.infile)
    if not args.invcf:
        logger.error("No vcf to filter specified!")
        raise ValueError("No vcf to filter specified!")
    # close the incoming handles to avoid double-opening the files
    #args.invcf.close()
    #args.outfile.close()
    io.extract_syntenic_from_vcf(df, args.invcf.name, args.outfile.name)


def plot(args):
    """Deprecated/internal, DO NOT USE


    """
    df = io.read_pff(args.infile)
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



