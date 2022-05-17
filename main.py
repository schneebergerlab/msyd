#!/usr/bin/python3
# in python, probably not worth cythonizing

import ingest
import pansyn

import logging
import logging.config

import argparse as ap
import os
import sys
import numpy as np

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""

if __name__ == "__main__": # testing
    import sys
    syris = []
    bams = []
    for fin in sys.argv[1:]:
        syris.append(fin + "syri.out")
        bams.append(fin + ".bam")

    df = pansyn.find_pansyn(syris, bams, sort=False)
    print(df)
    print("regions:", len(df))
    print("total lengths:", sum(map(lambda x: x[1][0].end-x[1][0].start,df.iterrows())))
    #df = graph_pansyn(sys.argv[1:], mode='overlap')
    #print(df)
    #print("regions:", len(df))
sys.exit(0)


parser = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
optional = parser._action_groups.pop()
required = parser.add_argument_group("Input Files")
required.add_argument("-c", dest="infile", help="File containing alignment coordinates", type=ap.FileType('r'), required=True)
#required.add_argument("-r", dest="ref", help="Genome A (which is considered as reference for the alignments). Required for local variation (large indels, CNVs) identification.", type=ap.FileType('r'))
#required.add_argument("-q", dest="qry", help="Genome B (which is considered as query for the alignments). Required for local variation (large indels, CNVs) identification.", type=ap.FileType('r'))
#required.add_argument("-d", dest="delta", help=".delta file from mummer. Required for short variation (SNPs/indels) identification when CIGAR string is not available", type=ap.FileType('r'))

other = parser.add_argument_group("Additional arguments")
other.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="T", choices=['T', 'S', 'B', 'P'])
other.add_argument('-f', dest='f', help='Filter out low quality alignments', default=True, action='store_false')
other.add_argument('-k', dest="keep", help="Keep intermediate output files", default=False, action="store_true")
other.add_argument('--dir', dest='dir', help="path to working directory (if not current directory). All files must be in this directory.", action='store')
other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="")
other.add_argument("--seed", dest="seed", help="seed for generating random numbers", type=int, default=1)
other.add_argument('--nc', dest="nCores", help="number of cores to use in parallel (max is number of chromosomes)", type=int, default=1)
other.add_argument('--novcf', dest="novcf", help="Do not combine all files into one output file", default=False, action="store_true")

# Parameters for identification of structural rearrangements
srargs = parser.add_argument_group("SR identification")
srargs.add_argument("--nosr", dest="nosr", help="Set to skip structural rearrangement identification", action="store_true", default=False)
srargs.add_argument("--invgaplen", dest="invgl", help="Maximum allowed gap-length between two alignments of a multi-alignment inversion.", type=int, default=500000)
srargs.add_argument("--tdgaplen", dest="tdgl", help="Maximum allowed gap-length between two alignments of a multi-alignment translocation or duplication (TD). Larger values increases TD identification sensitivity but also runtime.", type=int, default=500000)
srargs.add_argument("--tdmaxolp", dest="tdolp", help="Maximum allowed overlap between two translocations. Value should be in range (0,1].", type=float, default=0.8)
srargs.add_argument("-b", dest="bruteRunTime", help="Cutoff to restrict brute force methods to take too much time (in seconds). Smaller values would make algorithm faster, but could have marginal effects on accuracy. In general case, would not be required.", type=int, default=60)
srargs.add_argument("--unic", dest="TransUniCount", help="Number of uniques bps for selecting translocation. Smaller values would select smaller TLs better, but may increase time and decrease accuracy.", type=int, default=1000)
srargs.add_argument("--unip", dest="TransUniPercent", help="Percent of unique region requried to select translocation. Value should be in range (0,1]. Smaller values would allow selection of TDs which are more overlapped with \
 other regions.", type=float, default=0.5)
srargs.add_argument("--inc", dest="increaseBy", help="Minimum score increase required to add another alignment to translocation cluster solution", type=int, default=1000)
srargs.add_argument("--no-chrmatch", dest='chrmatch', help="Do not allow SyRI to automatically match chromosome ids between the two genomes if they are not equal", default=False, action='store_true')

# Parameters for identification of short variations
shvargs = parser.add_argument_group("ShV identification")
shvargs.add_argument("--nosv", dest="nosv", help="Set to skip structural variation identification", action="store_true", default=False)
shvargs.add_argument("--nosnp", dest="nosnp", help="Set to skip SNP/Indel (within alignment) identification", action="store_true", default=False)
# shvargs.add_argument("-align", dest="align", help="Alignment file to parse to show-snps for SNP/Indel identification", action="store_false", type=argparse.FileType("r"))
shvargs.add_argument("--all", help="Use duplications too for variant identification",  action="store_true", default=False)
shvargs.add_argument("--allow-offset", dest='offset', help='BPs allowed to overlap', default=5, type=int, action="store")
shvargs.add_argument('--cigar', dest="cigar", help="Find SNPs/indels using CIGAR string. Necessary for alignments generated using aligners other than nucmers", default=False, action='store_true')
shvargs.add_argument('-s', dest="sspath", help="path to show-snps from mummer", default="show-snps")
# shvargs.add_argument('--maxsvseq', dest="mss", help="Sets the upper limit on the size of SV for which bases would be outputted files. Currently, only affects insertions, deletions and HDRs. (-1 means no cut-off)", type=int, default=-1)

optional.add_argument("--lf", dest="log_fin", help="Name of log file", type=ap.FileType("w"), default="syri.log")
optional.add_argument("--log", dest="log", help="log level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARN"])

args = parser.parse_args()

logging.config.dictConfig({
'version': 1,
'disable_existing_loggers': False,
'formatters': {
    'log_file': {
        'format': "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s",
    },
    'stdout': {
        'format': "%(name)s - %(levelname)s - %(message)s",
    },
},
'handlers': {
    'stdout': {
        'class': 'logging.StreamHandler',
        'formatter': 'stdout',
        'level': 'WARNING',
    },
    'log_file': {
        'class': 'logging.FileHandler',
        'filename': args.log_fin.name,
        'mode': 'a',
        'formatter': 'log_file',
        'level': args.log,
    },
},
'loggers': {
    '': {
        'level': args.log,
        'handlers': ['stdout', 'log_file'],
    },
},
})


coords, chrlink = ingest.readCoords(args.infile.name, args.chrmatch, args.dir, args.prefix, args, args.cigar)
achrs = np.unique(coords.aChr).tolist()
bchrs = np.unique(coords.bChr).tolist()


print(coords)
