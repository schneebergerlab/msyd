#!/bin/python3

import pandas as pd
import pytest
import os

import pansyri.util as util

from pansyri.pansyn import *
from pansyri.classes.coords import Range, Pansyn
from pansyri.classes.cigar import Cigar


reffwd = set(['M', 'D', 'N', '=', 'X'])
qryfwd = set(['M', 'I', 'S', '=', 'X'])
cig_types = set(['M', '=', 'X', 'S', 'H', 'D', 'I', 'N'])
cig_aln_types = set(['M', 'X', '='])
cig_clips = set(['S', 'H', 'P', 'N']) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these

def test_pansyn_int():
    """
    Integration test testing the higher-order functionality of the pansyn module by validating the alignments.
    """
    ## init
    os.chdir('../../ampril_reduced/') # hardcoded for now to local test dataset
    syns, alns = util.parse_input_tsv('./full.tsv')
    
    # read in genome files
    genome_files = [syn.split('.')[0].split('_')[-1] + '.fasta.gz' for syn in syns]
    gens = {f.split('.')[0]:read_fasta(f) for f in genome_files} # key by org name
    refgen = read_fasta('./col.fasta.gz')
    
    # get pansyn df
    df = util.crosssyn_from_lists(syns, alns)

    ## do the validation
    for row in df:
        pan = row[1][0]
        refseq = refgen[pan.ref.chr][pan.ref.start:pan.ref.end]
        for org in pan.get_organisms():
            rng = pan.ranges_dict[org]
            cg = pan.cigars_dict[org]
            
            qryseq = gens[rng.org][rng.chr][rng.start:rng.end]

            progr = 0
            progq = 0
            for l, t in cg.pairs:
                assert(t in cig_types) # should be a valid CIGAR
                assert(t not in cig_clips)  # there shouldn't be any clipping

                refcmp = refseq[prog:progr+l]
                qrycmp = qryseq[prog:progq+l]
                
                # check concretely for matching types
                if t == '=':
                    assert(refcmp == qrycmp)
                elif t == 'X':
                    assert(refcmp != qrycmp)

                if t in reffwd:
                    progr += l
                if t in qryfwd:
                    progq += l

            # alignments should be total
            assert(progr == len(refseq))
            assert(progq == len(qryseq))
