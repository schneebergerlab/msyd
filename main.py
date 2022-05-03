#!/usr/bin/python3
# in python, probably not worth cythonizing

import ingest
import argparse as ap

import logging
import logging.config

import os
import sys

"""
This file serves as the main entrypoint for finding pansyntentic regions.
Experimental and WIP.
"""


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


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
optional = parser._action_groups.pop()
required = parser.add_argument_group("Input Files")
required.add_argument("-c", dest="infile", help="File containing alignment coordinates", type=argparse.FileType('r'), required=True)
#required.add_argument("-r", dest="ref", help="Genome A (which is considered as reference for the alignments). Required for local variation (large indels, CNVs) identification.", type=argparse.FileType('r'))
#required.add_argument("-q", dest="qry", help="Genome B (which is considered as query for the alignments). Required for local variation (large indels, CNVs) identification.", type=argparse.FileType('r'))
#required.add_argument("-d", dest="delta", help=".delta file from mummer. Required for short variation (SNPs/indels) identification when CIGAR string is not available", type=argparse.FileType('r'))

other = parser.add_argument_group("Additional arguments")
other.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="T", choices=['T', 'S', 'B', 'P'])
#other.add_argument('-f', dest='f', help='Filter out low quality alignments', default=True, action='store_false')
other.add_argument('-k', dest="keep", help="Keep intermediate output files", default=False, action="store_true")
other.add_argument('--dir', dest='dir', help="path to working directory (if not current directory). All files must be in this directory.", action='store')
other.add_argument("--prefix", dest="prefix", help="Prefix to add before the output file Names", type=str, default="")
other.add_argument("--seed", dest="seed", help="seed for generating random numbers", type=int, default=1)
other.add_argument('--nc', dest="nCores", help="number of cores to use in parallel (max is number of chromosomes)", type=int, default=1)
other.add_argument('--novcf', dest="novcf", help="Do not combine all files into one output file", default=False, action="store_true")

args = parser.parse_args()

coords, chrlink = readCoords(args.infile.name, args.chrmatch, args.dir, args.prefix, args, args.cigar)
achrs = np.unique(coords.aChr).tolist()
bchrs = np.unique(coords.bChr).tolist()
