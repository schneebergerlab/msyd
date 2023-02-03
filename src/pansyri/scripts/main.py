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

    parser = argparse.ArgumentParser(description="""
    Pansyn is a pansynteny identifier.
    Pansyn consists of a Cython library and a CLI interface.\n
    The CLI interface consists of multiple subcommands, described briefly below.\n
    For more information, see the documentation and subparser help messages accessed by calling pansyn [subparser] -h.
    """)
    parser.set_defaults(func=None, cores=1)

    subparsers = parser.add_subparsers()#description="See also pansyri [subparser] -h:") # title/description?
    # ordering parser
    order_parser = subparsers.add_parser("order",
        help="Determine a suitable ordering for plotting from a pansyn callset.",
        description="""
        Determine the optimal ordering of the supplied genomes for plotting using a clustering-based algorithm.
        The ordering is determined such that adjacent organisms share as many basepairs of pansynteny  as possible.
        """)
    order_parser.set_defaults(func=order)
    order_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")

    ## plotting subparser
    #plot_parser = subparsers.add_parser("plot", description="Prints a lengths df for plotting to stdout. Can be piped to a file and plotted with tests/plot.R .")
    #plot_parser.set_defaults(func=plot)
    #plot_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF or VCF file to read pansynteny information from.")

    # Pansyn calling argparser
    call_parser = subparsers.add_parser("call",
        help="Identify pansynteny from a set of alignments and syri calls to reference.",
        description="""
        Call Pansynteny in a set of genomes that have been aligned to a reference and processed with syri.\n
        Requires a tab-separated file listing for each organism the name that should be used, the path to the alignment and syri output files.\n
        Output can be saved either in Pansyri File Format (PFF) or VCF. VCF output does not preserve alignment information and cannot be used for some of the further processing!\n
        By default, Pansyn runs an exact pansynteny calling algorithm respecting the alignment information; for preliminary analyses, it might be suitable to use a faster, approximate algorithm.
        This can be done using some of the flags described below:
        """)
    call_parser.set_defaults(func=call)
    call_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="The .tsv file to read SyRI output and alignment files in from. For more details, see the Readme.")
    call_parser.add_argument("-o", dest='pff', type=argparse.FileType('wt'), help="Where to save the output in PFF format (see format.md). At least one of -o and -v must be specified!")
    call_parser.add_argument("-v", "--vcf", dest='vcf', type=argparse.FileType('wt'), help="Where to save the output in VCF format. At least one of -o and -v must be specified!")
    call_parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Pansyn cannot make effective use of more cores than the number of input organisms divided by two. Defaults to 4.", type=int, default=4)
    call_parser.add_argument("--core", dest='core', action='store_const', const=True, default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    call_parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL SyRI annotations. Yields more contiguous regions and faster runtime, but calls may not be exact to the base level.")
    call_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified.")
    call_parser.add_argument("-p", "--print", dest='print', action='store_true', default=False, help="print a subset of the output to stdout, for debugging.")

    ## Conversion parser
    #convert_parser = subparsers.add_parser("convert",
    #    help="Convert PFF to VCF files.",
    #    description="""
    #    Convert between different Pansynteny annotations file formats.
    #    Additional formats may be supported in the future.
    #    """)
    #convert_parser.set_defaults(func=convert)
    #convert_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    #convert_parser.add_argument("-t", dest='filetype', required=False, type=str, help="File format to convert to.")
    #convert_parser.add_argument("-o", dest='outfile', default='-', type=argparse.FileType('wt'), help="Where to store disk output. Defaults to stdout (specified with \"-\").")

    # Filter subparser
    filter_parser = subparsers.add_parser("filter",
        help="Filter a VCF file to only contain annotations in pansyntenic regions",
        description="""
        Used for filtering VCF files to only contain calls in pansyntenic regions.
        Can be run on pff files processed with pansyn view.
        """)
    filter_parser.add_argument("--vcf", dest='invcf', required=True, type=argparse.FileType('r'), help="The .vcf file to filter and write to -o.")
    filter_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    filter_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the filtered VCF.")


    # view subparser
    view_parser = subparsers.add_parser("view",
        help="Filter, convert or analyze existing PFF Files",
        description="""
        Used for filtering VCF files to only contain calls in pansyntenic regions for now.
        Additional functionality will be implemented later.
        """)
    view_parser.set_defaults(func=view)
    view_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    view_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the output. File format is determined automatically from the extension, but can be overridden by supplying any of the --o flags.")
    view_parser.add_argument("--opff", dest='filetype', action='store_const', const='pff', help="store output in PFF format")
    view_parser.add_argument("--opff-nocg", dest='filetype', action='store_const', const='pff-nocg', help="store output in PFF format, discarding cigar strings")
    view_parser.add_argument("--ovcf", dest='filetype', action='store_const', const='vcf', help="store output in VCF format, discarding cigar strings")

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

def view(args):
    logger.info(f"reading pansyn output from {args.infile.name}")
    df = io.read_pff(args.infile)
    if not args.filetype: # determine filetype if not present
        args.filetype = args.outfile.name.split(".")[-1]

    # do further processing here

    print(util.get_stats(df))
    # save
    logger.info(f"saving to {args.outfile.name} in {args.filetype} format")
    if args.filetype == 'pff':
        io.save_to_pff(df, args.outfile)
    elif args.filetype == 'vcf':
        io.save_to_vcf(df, args.outfile)
    elif args.filetype == 'pff-nocg' or args.filetype == 'pff-nocigar':
        io.save_to_pff(df, args.outfile, save_cigars=False)
    else:
        logger.error(f"Invalid filetype: {args.filetype}")

def filter(args):
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



