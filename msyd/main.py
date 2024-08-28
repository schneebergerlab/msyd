#!/usr/bin/env python3

import msyd # to import version
import msyd.util as util
import msyd.io as io
import msyd.vcf as vcf
import msyd.imputation as imputation
import msyd.multisyn as multisyn
import msyd.realignment as realignment
from msyd.coords import Range

import msyd.ordering as ordering

logger = util.CustomFormatter.getlogger(__name__)

import pandas as pd

import argparse
import sys
import os

"""
This file serves as the main entrypoint for the msyd CLI.
"""

def main():
    from msyd.util import logger
    from msyd import __version__ as msydv

    parser = argparse.ArgumentParser(description="""
    msyd is a tool for identifying and processing multisynteny.
    msyd consists of a Python library and a CLI interface.\n
    The CLI interface consists of multiple subcommands, described briefly below.\n
    For more information, see the documentation and subparser help messages accessed by calling msyd [subparser] -h.
    """)
    parser.set_defaults(func=None, cores=1)
    parser.add_argument('--version', action='version', version=msydv)

    subparsers = parser.add_subparsers()#description="See also msyd [subparser] -h:") # title/description?
    # ordering parser
    order_parser = subparsers.add_parser("order",
        help="Determine a suitable ordering for plotting from a multisynteny callset.",
        description="""
        Determine the optimal ordering of the supplied genomes for plotting using a clustering-based algorithm.
        The ordering is determined such that adjacent organisms share as many basepairs of multisynteny  as possible.
        """)
    order_parser.set_defaults(func=order)
    order_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PSF file to read multisynteny information from.")

    ## plotting subparser
    #plot_parser = subparsers.add_parser("plot", description="Prints a lengths df for plotting to stdout. Can be piped to a file and plotted with tests/plot.R .")
    #plot_parser.set_defaults(func=plot)
    #plot_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PSF or VCF file to read multisynteny information from.")

    ## Filter subparser
    #filter_parser = subparsers.add_parser("filter",
    #    help="Filter a VCF file to only contain annotations in multisyntenic regions",
    #    description="""
    #    Used for filtering VCF files to only contain calls in multisyntenic regions.
    #    Can be run on psf files processed with msyd view.
    #    """)
    #filter_parser.set_defaults(func=filter)
    #filter_parser.add_argument("--vcf", dest='invcf', required=True, type=argparse.FileType('r'), help="The .vcf file to filter and write to -o.")
    #filter_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PSF file to read multisynteny information from.")
    #filter_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the filtered VCF.")
    #filter_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="The reference to use for the synteny annotated in the output VCF")



    # Multisyn calling argparser
    call_parser = subparsers.add_parser("call",
        help="Identify multisynteny from a set of alignments and syri calls to reference.",
        description="""
        Call Multisynteny in a set of genomes that have been aligned to a reference and processed with syri.\n
        Requires a tab-separated file listing for each organism the name that should be used, the path to the alignment and syri output files.\n
        Output can be saved either in Population Synteny File Format (.psf) or VCF. VCF output does not preserve alignment information and cannot be used for some of the further processing!\n
        """)
    call_parser.set_defaults(func=call)
    call_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="The .tsv file to read SyRI output, alignment and VCF files in from. For more details, see the Readme.")
    call_parser.add_argument("-o", dest='psf', required=True, type=argparse.FileType('wt'), help="Where to save the output PSF file (see format.md)")
    call_parser.add_argument("-m", "--merge-vcf", dest='vcf', type=argparse.FileType('wt'), help="Merge the VCFs specified in the input table, store the merged VCF at the path specified. Does not currently work with --realign, as non-ref haplotypes do not have coordinates on the reference that VCF records can be fetched from.")
    call_parser.add_argument("-a", "--all", dest='all', action='store_true', default=False, help="Merge all VCF records instead of only records annotated in multisyntenic regions.")
    call_parser.add_argument("-x", "--complex", dest='no_complex', action='store_const', const=False, default=True, help="Do not filter the input VCFs to only contain SNPs and INDELs")
    call_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="Reference to use for the VCF output")
    call_parser.add_argument("--incremental", dest='incremental', type=argparse.FileType('r'), help="A PSF file containing a previous multisynteny callset to combine with the calls derived from the input TSV. Should contain CIGAR strings.")
    call_parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Multisyn cannot make effective use of more cores than the number of input organisms divided by two. Defaults to 1.", type=int, default=1)
    call_parser.add_argument("--core", dest='core', action='store_true', default=False, help="Call only core synteny. Improves runtime significantly, particularly on larger datasets.")
    call_parser.add_argument("--syn", "-s", dest='SYNAL', action='store_const', const=False, default=True, help="Use SYN instead of SYNAL SyRI annotations. Fast, but error-prone and inaccurate. Not recommended.")
    call_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .psf file. Has no effect when --syn is specified.")
    call_parser.add_argument("--realign", "-ali", dest='realign', action='store_true', default=False, help="After calling core and reference cross synteny, realign missing regions to identify non-reference synteny.")
    call_parser.add_argument("--pairwise", dest='pairwise', required=False, type=argparse.FileType('r'), help="Path to a TSV containing paths to full pairwise alignments that msyd will read in from disk during realignment if this parameter is passed. Otherwise, individual regions will be realigned on the fly with minimap2/mappy. This is useful if you already have pairwise alignments, or want to use a different aligner.")
    call_parser.add_argument("-p", "--print", dest='print', action='store_true', default=False, help="print a subset of the output to stdout, for debugging.")
    call_parser.add_argument("--impute", dest='impute', action='store_true', default=False, help="When processing small variants in a VCF, interpret the lack of a variant as identical to the reference genotype for that haplotype.")
    call_parser.add_argument("--workdir", "-w", dest='tmp', required=False, type=str, help="Path to a working directory to be used for storing temporary files. If the path does not exist, it will be created!")
    call_parser.add_argument("--min-realign", dest="min_realign", help="Minimum region size to realign, in bp. Default 150 bp.", type=int, default=-1)
    call_parser.add_argument("--min-syn-id", dest="min_syn_id", help="% Identity required for a region to be called as syntenic during the realignment step. Default 80%.", type=int, default=80)
    call_parser.add_argument("--max-realign", dest="max_realign", help="Maximum number of realignment steps to perform. Default 0 (unlimited).", type=int, default=-1)
    call_parser.add_argument("--minimap-preset", dest="mp_preset", help="minimap2 alignment preset to use. Default 'asm20'.", type=str, default="asm20")

    # view subparser
    view_parser = subparsers.add_parser("view",
        help="Filter, convert or analyze existing PSF Files",
        description="""
        Used for filtering VCF files to only contain calls in multisyntenic regions for now.
        Additional functionality will be implemented later.
        """)
    view_parser.set_defaults(func=view)
    view_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PSF file to read multisynteny information from.")
    view_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to store the output. File format is determined automatically from the extension, but can be overridden by supplying any of the --o flags.")
    view_parser.add_argument("-e", dest='expr', action='store', type=str, help="Expression to use for filtering the multisyntenic regions. This is done before --intersect is evaluated if also supplied")
    view_parser.add_argument("-p", dest='print', action='store_const', const=10, help="Print the first 10 regions after filtering, mainly for debugging")
    view_parser.add_argument("-r", "--reference", dest='ref', type=argparse.FileType('r'), help="If saving to VCF, the reference to use can be specified with this flag")
    view_parser.add_argument("--intersect", dest='intersect', type=argparse.FileType('r'), help="VCF File to intersect with the PSF file given with -i. Will only keep annotations within multisyntenic regions")
    view_parser.add_argument("--impute", dest='impute', action='store_true', default=False, help="When processing small variants in a VCF, interpret the lack of a variant as identical to the reference genotype for that haplotype.")

    view_parser.add_argument("--opsf", dest='filetype', action='store_const', const='psf', help="store output in PSF format")
    view_parser.add_argument("--opsf-nocg", dest='filetype', action='store_const', const='psf-nocg', help="store output in PSF format, discarding cigar strings")
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
        help="Iteratively realign a set of genomes based on a PSF file",
        description="""
        Exposes the realignment functionality in msyd call directly.
        Useful for realigning only a specific region by prefiltering the PSF.
        """)
    realign_parser.set_defaults(func=realign)
    realign_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PSF file to read multisynteny information from.")
    realign_parser.add_argument("-o", dest='outfile', required=True, type=argparse.FileType('wt'), help="Where to save the output PSF file (see format.md)")
    realign_parser.add_argument("-t", dest='tsvfile', required=True, type=argparse.FileType('r'), help="TSV containing the sample names and path to genome fastas.")
    realign_parser.add_argument("-p", "--pairwise", dest='pairwise', required=False, type=argparse.FileType('r'), help="Path to a TSV containing paths to full pairwise alignments that msyd will read in from disk if this parameter is passed. Otherwise, individual regions will be realigned on the fly with minimap2/mappy. This is useful if you already have pairwise alignments, or want to use a different aligner.")
    realign_parser.add_argument("--workdir", "-w", dest='tmp', required=False, type=str, help="Path to a working directory to be used for storing temporary files. If the path does not exist, it will be created!")
    realign_parser.add_argument("--no-cigars", dest='cigars', action='store_const', const=False, default=True, help="Don't store CIGAR strings in the saved .psf file. Has no effect when --syn is specified.")
    realign_parser.add_argument("--min-realign", dest="min_realign", help="Minimum region size to realign, in bp. Default 100 bp.", type=int, default=-1)
    realign_parser.add_argument("--min-syn-id", dest="min_syn_id", help="% Identity required for a region to be called as syntenic during the realignment step. Default 80%.", type=int, default=80)
    realign_parser.add_argument("--max-realign", dest="max_realign", help="Maximum number of realignment steps to perform. Default 0 (unlimited).", type=int, default=-1)
    realign_parser.add_argument("--minimap-preset", dest="mp_preset", help="minimap2 alignment preset to use. Default 'asm20'.", type=str, default="asm20")


    stats_parser = subparsers.add_parser("stats",
        help="Compute some statistics on a PSF file",
        description="""
        Computes some basic statistics on a PSF file.
        Useful as input for plotting or to get a feel for the dataset.
        """)
    stats_parser.set_defaults(func=stats)
    stats_parser.add_argument("-i", dest='infile', required=True, type=argparse.FileType('r'), help="PSF file to read multisynteny information from.")
    stats_parser.add_argument("-o", dest='outfile', default='-', type=argparse.FileType('wt'), help="Where to send the statistics to. Default stdout.")
    stats_parser.add_argument("--separator", "-s", dest="sep", help="Separator to use for printing the stats. Default is tab (for TSV), set to ',' for CSV.", type=str, default="\t")
    stats_parser.add_argument("-p", "--prefix", dest='siprefix', action='store_true', default=False, help="Whether to attach SI prefixes to the output for human readability. If not supplied, print exact numbers.")
    stats_parser.add_argument("-a", "--aggregate", dest='agg', action='store_true', default=False, help="If passed, will report summary statistics for all haplotypes instead of by organism.")
    stats_parser.add_argument("--no-header", dest='header', action='store_false', default=True, help="If passed, msyd will not print a header for the CSV.")
    #stats_parser.add_argument("-r", "--reference", dest='agg', action='store_true', default=False, help="If passed, will report summary statistics")

    #fact_parser = subparsers.add_parser("fact",
    #    help="Give a fact about birds or non-birds!",
    #    description="""
    #        Written in a duck-typed language!
    #    """)
    order_parser.set_defaults(func=order)

    args = parser.parse_args()
    if args.func:
        logger.info("Starting msyd.")
        args.func(args)
        logger.info("Finished running msyd. Have a nice day!")
    else:
        logger.info("No subcommand specified, printing help message.")
        parser.print_help()
    return
# END

def call(args):
    # import msyd  # to import version
    # from msyd.script.io import find_multisyn
    import msyd.io as io
    import msyd.imputation as imputation
    import msyd.realignment as realignment
    import msyd.intersection as intersection
    from msyd.coords import Range
    import msyd.ordering as ordering

    import msyd.util as util
    logger = util.CustomFormatter.getlogger("call")

    logger.info("Starting msyd call")

    qrynames, syns, alns, vcfs, fastas = util.parse_input_tsv(args.infile)
    # find reference synteny
    #syndicts = intersection.find_multisyn(qrynames, syns, alns, only_core=args.core, SYNAL=args.SYNAL, base=args.incremental)
    syndict = intersection.prepare_input(qrynames, syns, alns, cores=args.cores, SYNAL=args.SYNAL, base=args.incremental)
    logger.info("Read input files")

    syndict = intersection.process_syndicts(syndict, cores=args.cores)
    logger.info("Intersected synteny")

    if args.realign:
        # read in full pairwise alns if supplied
        if args.pairwise:
            alndict = io.read_alnsfile(args.pairwise)
        #TODO directly do in call to realignment

        # use reference synteny as base to identify all haplotypes
        df = realignment.realign(df, qrynames, fastas, MIN_REALIGN_LEN=args.min_realign, MIN_SYN_ID=args.min_syn_id, MAX_REALIGN=args.max_realign, mp_preset=args.mp_preset, ncores=args.cores, pairwise=alndict if args.pairwise else None)

        # garb = realign(df, qrynames, fastas, MIN_REALIGN_LEN=args.min_realign, MAX_REALIGN=args.max_realign, mp_preset=args.mp_preset, ncores=args.cores, cwd=args.tmp)
        # realign(syns, qrynames, fastas, MIN_REALIGN_LEN=None, MAX_REALIGN=None, mp_preset='asm5'):

    if args.tmp:
        if not os.path.isdir(args.tmp):
            os.makedirs(args.tmp)
        TMPDIR = args.tmp

    if args.print:
        logger.info("Printing sample head to STDOUT")
        print(df.head())

    print(util.get_map_stats(df))

    # save output
    logger.info(f"Saving msyd calls to PSF at {args.psf.name}")
    io.save_to_psf(df, args.psf, save_cigars=args.cigars)

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
            logger.info("Pre-filtering VCFs to multisyntenic regions")
            vcfs = vcf.filter_vcfs(df, vcfs, ref, no_complex=args.no_complex, add_syn_anns=False, impute_ref=args.impute)
            

        logger.info(f"Filtered files: {vcfs}")

        tmpfile = util.gettmpfile()
        logger.info(f"Merging VCFs, saving to {tmpfile}")
        vcf.reduce_vcfs(vcfs, tmpfile)

        if args.impute:
            logger.info(f"Imputing reference genotypes in syntenic regions, saving to {args.vcf.name}")
            vcf.extract_syntenic_from_vcf(df, tmpfile, args.vcf.name, no_complex=args.no_complex, add_syn_anns=True, impute_ref=args.impute)
        else:
            logger.info(f"Adding multisynteny annotations, saving to {args.vcf.name}")
            vcf.add_syn_anns_to_vcf(df, tmpfile, args.vcf.name, ref=ref) 

    logger.info(f"Finished running msyd call, output saved to {args.psf.name}.")

def merge(args):
    import msyd.io as io

    import msyd.util as util
    logger = util.CustomFormatter.getlogger("merge")

    # temporary function to better test the vcf merging functionality
    logger.info(f"Merging {args.vcfs} to {args.outfile.name}")
    vcf.reduce_vcfs(args.vcfs, args.outfile.name)
    logger.info(f"Finished running msyd merge, output saved to {args.outfile.name}.")

# call the plotsr ordering functionality on a set of organisms described in the .tsv
def order(args):
    import msyd.io as io
    import msyd.ordering as ordering

    import msyd.util as util
    logger = util.CustomFormatter.getlogger("order")

    df = io.read_psf(args.infile)
    print(ordering.order_hierarchical(df, orgs=None, score_fn=ordering.syn_score))
    logger.info("Finished running msyd order")

def view(args):
    import msyd.io as io
    import msyd.util as util
    logger = util.CustomFormatter.getlogger("view")
    logger.info(f"reading multisynteny from {args.infile.name}")
    df = io.read_psf(args.infile)
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
        vcf.extract_syntenic_from_vcf(df, args.intersect.name, args.outfile.name, ref=args.ref.name if args.ref else None, impute_ref=args.impute)
        return # has been saved already

    # save
    logger.info(f"Writing to {args.outfile.name} in {args.filetype} format")
    if args.filetype == 'psf':
        io.save_to_psf(df, args.outfile)
    elif args.filetype == 'vcf':
        io.save_to_vcf(df, args.outfile, args.ref.name if args.ref else None)
    elif args.filetype == 'psf-nocg' or args.filetype == 'psf-nocigar':
        io.save_to_psf(df, args.outfile, save_cigars=False)
    else:
        logger.error(f"Couldn't determine filetype for {args.filetype}")
        return
    logger.info(f"Finished running msyd view, output saved to {args.outfile.name}.")


def realign(args):
    import msyd.io as io
    import msyd.realignment as realignment

    import msyd.util as util
    logger = util.CustomFormatter.getlogger("realign")

    # read in full pairwise alns if supplied
    if args.pairwise:
        alndict = io.read_alnsfile(args.pairwise)

    logger.info(f"realigning from {args.infile.name}, taking genome files from {args.tsvfile.name}")
    qrynames, syris, alns, vcfs, fastas = util.parse_input_tsv(args.tsvfile)
    syns = io.read_psf(args.infile)
    resyns = realignment.realign(syns, qrynames, fastas, MIN_REALIGN_LEN=args.min_realign, MIN_SYN_ID=args.min_syn_id, MAX_REALIGN=args.max_realign, pairwise=alndict if args.pairwise else None)
    print(util.get_stats(resyns))

    logger.info(f"Saving to {args.outfile.name} in PSF format.")
    io.save_to_psf(resyns, args.outfile, save_cigars=args.cigars)
    logger.info(f"Finished running msyd realign, output saved to {args.outfile.name}.")

def stats(args):
    import msyd.io as io
    import msyd.util as util
    logger = util.CustomFormatter.getlogger("stats")

    logger.info(f"Reading from {args.infile.name}.")
    syns = io.read_psf(args.infile)
    #print(util.get_stats(resyns), file=args.outfile)
    if args.agg:
        print(util.get_stats(syns), file=args.outfile)
    else:
        print(util.lensdict_to_table(util.tabularize_lens_byorg(syns), sep=args.sep, si=args.siprefix, header=args.header), file=args.outfile)
    logger.info(f"Finished running msyd stats, output printed to {args.outfile.name}.")

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
    import msyd.io as io
    df = io.read_psf(args.infile)
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



