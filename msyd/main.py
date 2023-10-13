#!/usr/bin/python3
# in python, probably not worth cythonizing

import msyd # to import version
import msyd.util as util
import msyd.io as io
import msyd.imputation as imputation
import msyd.pansyn as pansyn
import msyd.realignment as realignment
from msyd.classes.coords import Range

import msyd.scripts.ordering as ordering

logger = util.CustomFormatter.getlogger(__name__)

import pandas as pd

import argparse
import sys
import os

"""
This file serves as the main entrypoint for the msyd CLI.
"""

def main(argv):

    parser = argparse.ArgumentParser(description="""
    msyd is a tool for identifying and processing pansynteny.
    msyd consists of a Python library and a CLI interface.\n
    The CLI interface consists of multiple subcommands, described briefly below.\n
    For more information, see the documentation and subparser help messages accessed by calling msyd [subparser] -h.
    """)
    parser.set_defaults(func=None, cores=1)
    parser.add_argument('--version', action='version', version=msyd.__version__)

    subparsers = parser.add_subparsers()#description="See also msyd [subparser] -h:") # title/description?
    # ordering parser
    order_parser = subparsers.add_parser("order",
        help="Determine a suitable ordering for plotting from a pansynteny callset.",
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
    #    Can be run on pff files processed with msyd view.
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
    call_parser.add_argument("-m", "--merge-vcf", dest='vcf', type=argparse.FileType('wt'), help="Merge the VCFs specified in the input table, store the merged VCF at the path specified.")
    call_parser.add_argument("-a", "--all", dest='all', action='store_const', const=True, default=False, help="Merge all VCF records instead of only records annotated in pansyntenic regions.")
    call_parser.add_argument("-x", "--complex", dest='no_complex', action='store_const', const=False, default=True, help="Do not filter the input VCFs to only contain SNPs and INDELs")
    call_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="Reference to use for the VCF output")
    call_parser.add_argument("--incremental", dest='incremental', type=argparse.FileType('r'), help="A PFF file containing a previous pansynteny callset to combine with the calls derived from the input TSV. Should contain CIGAR strings.")
    call_parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Pansyn cannot make effective use of more cores than the number of input organisms divided by two. Defaults to 1.", type=int, default=1)
    call_parser.add_argument("--core", dest='core', action='store_const', const=True, default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    call_parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL SyRI annotations. Yields more contiguous regions and faster runtime, but calls may not be exact to the base level.")
    call_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified.")
    call_parser.add_argument("--realign", "-ali", dest='realign', action='store_const', const=True, default=False, help="After calling core and reference cross synteny, realign missing regions to identify non-reference synteny.")
    call_parser.add_argument("-p", "--print", dest='print', action='store_true', default=False, help="print a subset of the output to stdout, for debugging.")
    call_parser.add_argument("--workdir", "-w", dest='tmp', required=False, type=str, help="Path to a working directory to be used for storing temporary files. If the path does not exist, it will be created!")
    call_parser.add_argument("--min-realign", dest="min_realign", help="Minimum region size to realign, in bp. Default 150 bp.", type=int, default=-1)
    call_parser.add_argument("--max-realign", dest="max_realign", help="Maximum number of realignment steps to perform. Default 0 (unlimited).", type=int, default=-1)

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

    merge_parser = subparsers.add_parser("merge",
        help="Merge different VCFs",
        description="""
        Exposes the optional VCF merging functionality in msyd call directly.
        Mainly for testing and debugging purposes
        """)
    merge_parser.set_defaults(func=merge)
    merge_parser.add_argument("-v", dest='vcfs', nargs='+', required=True, type=argparse.FileType('r'), help="The VCF files to merge.")
    merge_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the merged VCF.")

    realign_parser = subparsers.add_parser("realign",
        help="Iteratively realign a set of genomes based on a PFF file",
        description="""
        Exposes the realignment functionality in msyd call directly.
        Useful for realigning only a specific region by prefiltering the PFF.
        """)
    realign_parser.set_defaults(func=realign)
    realign_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PFF file to read pansynteny information from.")
    realign_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to save the output PFF file (see format.md)")
    realign_parser.add_argument("-t", dest='tsvfile', required=True, type=argparse.FileType('r'), help="TSV containing the sample names and path to genome fastas.")
    realign_parser.add_argument("--workdir", "-w", dest='tmp', required=False, type=str, help="Path to a working directory to be used for storing temporary files. If the path does not exist, it will be created!")
    realign_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .pff file. Has no effect when --syn is specified.")
    realign_parser.add_argument("--min-realign", dest="min_realign", help="Minimum region size to realign, in bp. Default 150 bp.", type=int, default=-1)
    realign_parser.add_argument("--max-realign", dest="max_realign", help="Maximum number of realignment steps to perform. Default 0 (unlimited).", type=int, default=-1)

    fact_parser = subparsers.add_parser("fact",
        help="Give a fact about birds or non-birds!",
        description="""
            Written in a duck-typed language!
        """)
    order_parser.set_defaults(func=order)

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)
        logger.info("Finished running msyd. Have a nice day!")
    else:
        logger.info("No subcommand specified, printing help message.")
        parser.print_help()

def call(args):
    qrynames, syns, alns, vcfs, fastas = util.parse_input_tsv(args.infile)
    # find reference synteny
    df = pansyn.find_multisyn(qrynames, syns, alns, only_core=args.core, SYNAL=args.SYNAL, base=args.incremental)
    if args.realign:
        # use reference synteny as base to identify all haplotypes
        df = realignment.realign(df, qrynames, fastas, MIN_REALIGN_THRESH=args.min_realign, MAX_REALIGN=args.max_realign)

    if args.tmp:
        if not os.path.isdir(args.tmp):
            os.makedirs(args.tmp)
        util.TMPDIR = args.tmp

    if args.print:
        logger.info("Printing sample head to STDOUT")
        print(df.head())

    print(util.get_stats(df))

    # save output
    logger.info(f"Saving msyd calls to PFF at {args.pff.name}")
    io.save_to_pff(df, args.pff, save_cigars=args.cigars)

    # if specified, merge the VCFs
    if args.vcf:
        logger.info(f"Merging VCFs into {args.vcf.name}")
        ref = None
        if args.ref:
            logger.info("Reading in Reference")
            ref = io.readfasta(args.ref.name)
        else:
            logger.warning("No reference specified. Specifying a reference is highly recommended to obtain standards-conforming VCFs!")

        if not args.all:
            logger.info("Pre-filtering VCFs to pansyntenic regions")
            vcfs = io.filter_vcfs(df, vcfs, ref, no_complex=args.no_complex, add_syn_anns=False)

        logger.info(vcfs)

        tmpfile = util.gettmpfile()
        logger.info(f"Merging VCFs, saving to {tmpfile}")
        io.reduce_vcfs(vcfs, tmpfile)

        logger.info(f"Adding pansynteny annotations, saving to {args.vcf.name}")
        io.add_syn_anns_to_vcf(df, tmpfile, args.vcf.name, ref=ref) 

    logger.info(f"Finished running msyd call, output saved to {args.pff.name}.")

def merge(args):
    # temporary function to better test the vcf merging functionality
    logger.info(f"Merging {args.vcfs} to {args.outfile.name}")
    io.reduce_vcfs(args.vcfs, args.outfile.name)
    logger.info(f"Finished running msyd merge, output saved to {args.outfile.name}.")

# call the plotsr ordering functionality on a set of organisms described in the .tsv
def order(args):
    df = io.read_pff(args.infile)
    print(ordering.order_hierarchical(df, orgs=None, score_fn=ordering.syn_score))
    logger.info("Finished running msyd order")

def view(args):
    logger.info(f"reading pansynteny from {args.infile.name}")
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
        return
    logger.info(f"Finished running msyd view, output saved to {args.outfile.name}.")


def realign(args):
    logger.info(f"realigning from {args.infile.name}, taking genome files from {args.tsvfile.name}")
    qrynames, syris, alns, vcfs, fastas = util.parse_input_tsv(args.tsvfile)
    syns = io.read_pff(args.infile)
    resyns = realignment.realign(syns, qrynames, fastas, MIN_REALIGN_THRESH=args.min_realign, MAX_REALIGN=args.max_realign)
    print(util.get_stats(resyns))
    logger.info(f"Saving to {args.outfile.name} in PFF format.")
    io.save_to_pff(resyns, args.outfile, save_cigars=args.cigars)
    logger.info(f"Finished running msyd realign, output saved to {args.outfile.name}.")

def fact(args):
    # print out a fact!
    import random
    with open("PATHTOJOKEFILE", 'r') as f:
        facts = []
        seq = ''
        #TODO implement tagging
        for line in f.readlines():
            if line[0] == '#':
                continue
            elif line == '===':
                facts.append(seq)
                seq = ''
        print(random.choice(facts))

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



