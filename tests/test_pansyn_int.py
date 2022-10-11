#!/bin/python3

import pandas as pd
import pytest
import os
import gzip
import platform

import pansyri.util as util

from pansyri.pansyn import *
from pansyri.classes.coords import Range, Pansyn
from pansyri.classes.cigar import Cigar


reffwd = set(['M', 'D', 'N', '=', 'X'])
qryfwd = set(['M', 'I', 'S', '=', 'X'])
cig_types = set(['M', '=', 'X', 'S', 'H', 'D', 'I', 'N'])
cig_aln_types = set(['M', 'X', '='])
cig_clips = set(['S', 'H', 'P', 'N']) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these


def read_fasta(f):
    """Helper function to read a fasta file given as a string path into a dictionary of chromosomes.
    """
    ret = {}
    with gzip.open(f, 'rt') as fin:
        key = ''
        strbuf = ''
        for line in fin:
            if line[0] == '>':
                if key:
                    ret[key] = strbuf
                    strbuf = ''
                key = line[1:].strip()
            else:
                strbuf += line.strip()
        ret[key] = strbuf # add last chr

    return ret


def test_pansyn_int():
    """
    Integration test testing the higher-order functionality of the pansyn module by validating the alignments.
    """
    ## init
    # on the cluster, go into full ampril, locally go into ampril_reduced
    os.chdir('../../ampril_reduced' if platform.node() == 'matmobile' else '../../data/ampril/')
    syns, alns = util.parse_input_tsv_path('./full.tsv')
    
    # read in genome files
    genome_files = [aln.split('.')[0].split('_')[-1] + '.filtered.fa.gz' for aln in alns]
    gens = {f.split('.')[0]:read_fasta(f) for f in genome_files} # key by org name
    refgen = read_fasta('./col.filtered.fa.gz')
    cnt = 0
    totlen = 0
    
    # get pansyn df
    df = util.crosssyn_from_lists(syns, alns)

    ## do the validation
    for row in df.iterrows():
        rowcnt = 0
        pan = row[1][0]
        refseq = refgen[pan.ref.chr][pan.ref.start:pan.ref.end +1]
        #print(pan)
        for org in pan.get_organisms():
            rng = pan.ranges_dict[org]
            cg = pan.cigars_dict[org]
            #print(org, rng, cg)
            
            qryseq = gens[rng.org][rng.chr][rng.start:rng.end +1]

            progr = 0
            progq = 0
            for l, t in cg.pairs:
                assert(t in cig_types) # should be a valid CIGAR
                assert(t not in cig_clips)  # there shouldn't be any clipping

                refcmp = refseq[progr:progr+l]
                qrycmp = qryseq[progq:progq+l]

                #print(progr, progq, l, t)
                
                # check concretely for matching types
                if t == '=':
                    rowcnt += sum(map(lambda x: 1 if x[0] != x[1] else 0, zip(refcmp, qrycmp)))
                    #assert(refcmp == qrycmp)
                elif t == 'X':
                    rowcnt += sum(map(lambda x: 1 if x[0] == x[1] else 0, zip(refcmp, qrycmp)))
                    #assert(refcmp != qrycmp)

                if t in reffwd:
                    progr += l
                if t in qryfwd:
                    progq += l

            # alignments should be total
            assert(progr == len(refseq))
            assert(progq == len(qryseq))

        cnt += rowcnt
        totlen += len(pan.ref)
        #print(cnt/len(pan.ref), cnt, len(pan.ref), pan)
    print(cnt/totlen, cnt, totlen)
    raise ValueError()
