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

    ## Filter subparser
    #filter_parser = subparsers.add_parser("filter",
    #    help="Filter a VCF file to only contain annotations in pansyntenic regions",
    #    description="""
    #    Used for filtering VCF files to only contain calls in pansyntenic regions.
    #    Can be run on pff files processed with pansyn view.
    #    """)
    #filter_parser.set_defaults(func=filter)
    #filter_parser.add_argument("--vcf", dest='invcf', required=True, type=argparse.FileType('r'), help="The .vcf file to filter and write to -o.")
    #filter_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    #filter_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the filtered VCF.")
    #filter_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="The reference to use for the synteny annotated in the output VCF")



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
    call_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="The .tsv file to read SyRI output, alignment and VCF files in from. For more details, see the Readme.")
    call_parser.add_argument("-o", dest='pff', required=True, type=argparse.FileType('wt'), help="Where to save the output PFF file (see format.md)")
    call_parser.add_argument("-v", "--vcf", dest='vcf', type=argparse.FileType('wt'), help="Where to save the merged VCF.")
    call_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="Reference to use for the VCF output")
    call_parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Pansyn cannot make effective use of more cores than the number of input organisms divided by two. Defaults to 4.", type=int, default=4)
    call_parser.add_argument("--core", dest='core', action='store_const', const=True, default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    call_parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL SyRI annotations. Yields more contiguous regions and faster runtime, but calls may not be exact to the base level.")
    call_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified.")
    call_parser.add_argument("-p", "--print", dest='print', action='store_true', default=False, help="print a subset of the output to stdout, for debugging.")

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
    view_parser.add_argument("-e", dest='expr', action='store', type=str, help="Expression to use for filtering the pansyntenic regions. This is done before --intersect is evaluated if also supplied")
    view_parser.add_argument("-p", dest='print', action='store_const', const=10, help="Print the first 10 regions after filtering, mainly for debugging")
    view_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="If saving to VCF, the reference to use can be specified with this flag")
    view_parser.add_argument("--intersect", dest='intersect', type=argparse.FileType('r'), help="VCF File to intersect with the PFF file given with -i. Will only keep annotations within pansyntenic regions")

    view_parser.add_argument("--opff", dest='filetype', action='store_const', const='pff', help="store output in PFF format")
    view_parser.add_argument("--opff-nocg", dest='filetype', action='store_const', const='pff-nocg', help="store output in PFF format, discarding cigar strings")
    view_parser.add_argument("--ovcf", dest='filetype', action='store_const', const='vcf', help="store output in VCF format, discarding cigar strings")

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)
    else:
        logger.info("No subcommand specified, priting help message.")
        parser.print_help()

def call(args):
    qrynames, syns, alns, vcfs = util.parse_input_tsv(args.infile)
    df = pansyn.find_multisyn(qrynames, syns, alns, only_core=args.core, SYNAL=args.SYNAL)
    if args.print:
        logger.info("Printing sample head to STDOUT")
        print(df.head())

    print(util.get_stats(df))

    # if specified, merge the VCFs
    if args.vcf:
        logger.info("Pre-filtering VCFs to syntenic regions")
        vcfs_filtered = io.prefilter(vcfs)
        logger.info("Merging VCFs, saving to {args.vcf.name}")
        io.reduce_vcfs(vcfs_filtered, args.vcf.name)

    # save output
    logger.info("Saving pansyn calls to PFF at {args.pff.name}")
    io.save_to_pff(df, args.pff, save_cigars=args.cigars)
    #if args.vcf:
    #    io.save_to_vcf(df, args.vcf, args.ref.name if args.ref else None, cores=args.cores)

# call the plotsr ordering functionality on a set of organisms described in the .tsv
def order(args):
    df = io.read_pff(args.infile)
    print(ordering.order_hierarchical(df, orgs=None, score_fn=ordering.syn_score))

def view(args):
    logger.info(f"reading pansyn output from {args.infile.name}")
    df = io.read_pff(args.infile)
    if not args.filetype: # determine filetype if not present
        args.filetype = args.outfile.name.split(".")[-1]
        logger.info(f"No output filetype specified - guessing from OUTFILE extension")

    # do further processing here
    if args.expr:
        logger.info(f"Applying filter: {args.expr}")
        df = util.apply_filtering(df, args.expr)

    if args.print:
        print(df.head(args.print))
    print(util.get_stats(df))

    if args.intersect:
        logger.info(f"Writing intersection to {args.outfile.name} as VCF")
        io.extract_syntenic_from_vcf(df, args.intersect.name, args.outfile.name, ref=args.ref.name if args.ref else None)
        return # has been saved already

    # save
    logger.info(f"Writing to {args.outfile.name} in {args.filetype} format")
    if args.filetype == 'pff':
        io.save_to_pff(df, args.outfile)
    elif args.filetype == 'vcf':
        io.save_to_vcf(df, args.outfile, args.ref.name if args.ref else None)
    elif args.filetype == 'pff-nocg' or args.filetype == 'pff-nocigar':
        io.save_to_pff(df, args.outfile, save_cigars=False)
    else:
        logger.error(f"Invalid filetype: {args.filetype}")

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



