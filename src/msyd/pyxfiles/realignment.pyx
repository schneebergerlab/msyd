#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3


import pandas as pd
import numpy as np
import mappy as mp
import pysam
import intervaltree

from collections import deque
import os
import multiprocessing
from functools import partial

from syri.synsearchFunctions import syri, mergeOutputFiles, outSyn
from syri.tdfunc import getCTX
from syri.writeout import getsrtable

import msyd.util as util


cdef int MIN_REALIGN_THRESH = 100
logger = util.CustomFormatter.getlogger(__name__)

cpdef realign(syns, qrynames, fastas):
    return process_gaps(syns, qrynames, fastas)

cdef process_gaps(syns, qrynames, fastas):
    """
    Function to find gaps between two coresyn regions and realign them to a new reference.
    Discovers all crosssynteny, hopefully.
    
    :arguments: A DataFrame with core and crosssyn regions called by find_pansyn and the sample genomes
    :returns: A DataFrame with the added non-reference crosssynteny
    """
    # init stuff
    ret = pd.DataFrame()
    n = len(qrynames)
    if not n == len(fastas):
        logger.error(f"More/less query names than fastas passed to process_gaps: {qrynames}, {fastas}!")
        raise ValueError("Wrong number of fastas!")
    
    # load fasta files
    fastas = {qrynames[i]: pysam.FastaFile(fastas[i]) for i in range(len(qrynames))}

    # iterate through each gap between coresyn blocks
    # call the alignment/ functionality and merge   

    syniter = syns.iterrows()
    try:
        syn = next(syniter)[1][0]

        # skip to first core
        while syn.get_degree() < n:
            syn = next(syniter)[1][0]
        old = syn

        while True:
            # find block between two coresyn regions
            crosssyns = []

            while syn.get_degree() < n:
                crosssyns.append(syn)
                syn = next(syniter)[1][0]
            # syn must be core now

            # start and end of the non-ref region, on the reference
            end = syn.ref.start
            start = old.ref.end

            # preemptively skip regions too small on the reference, if present
            if end - start < MIN_REALIGN_THRESH:
                syn = next(syniter)[1][0]
                continue

            # between chromosomes, there isn't a single gap
            if syn.ref.chr != old.ref.chr:
                old = syn
                syn = next(syniter)[1][0]
                continue


            # Block has been extracted and is long enough;
            # extract appropriate sequences, respecting crossyn

            #print(old, syn, syn.ref.start - old.ref.end)
            #print(crosssyns)

            # construct a dictionary containing for each sample a list of intervals that should be realigned
            mappingtrees = dict()
            seqdict = dict()
            for org in qrynames:
                chr = syn.ranges_dict[org].chr # chr must always stay the same
                pos = 0
                offset = old.ranges_dict[org].end # offset of the index in the new sequence to the old genome
                tree = intervaltree.IntervalTree()
                seq = ''
                fasta = fastas[org]
                for crosssyn in crosssyns:
                    if not org in crosssyn.ranges_dict: # no need to add this
                        continue
                    synrng = crosssyn.ranges_dict[org]
                    l = synrng.start - offset # len of the region to be added
                    if l < MIN_REALIGN_THRESH: # the allowed region is too small to add to realign
                        offset = synrng.end # skip till the end
                        continue

                    #print(l, offset, seq)
                    # add to the intervaltree
                    tree[pos:pos+l] = offset - pos # subtract starting position of the interval in collated sequence, as it will be readded later
                    seq += fasta.fetch(region=chr, start=offset, end=offset+l)
                    offset += l
                    pos += l

                # see if theres any sequence left to realign after processing the crosssyn regions
                l = syn.ranges_dict[org].start - offset
                if l >= MIN_REALIGN_THRESH:
                    tree[pos:pos+l] = offset
                    seq += fasta.fetch(region=chr, start=offset, end=offset+l)

                if tree and seq:
                    mappingtrees[org] = tree
                    seqdict[org] = seq
                elif seq or tree:
                    logger.error(f"Non-empty Tree with empty seq or the other way round: {tree}, {seq}")
                else:
                    pass
                    #logger.info(f"Leaving out {org}")

            if not seqdict: # if all sequences have been discarded, skip realignment
                logger.info("Not aligning, not enough non-reference sequence found!")
                old = syn
                syn = next(syniter)[1][0]
                continue

            # choose a reference as the sample containing the most non-crosssynteny
            ref = max(map(lambda x: (len(x[1]), x[0]), seqdict.items()))[1]
            print('ref:', ref)
            print('On ref:', syn.ref.chr, start, end, end - start)
            print({org:len(seq) for org, seq in seqdict.items()})
            #for org in seqdict:
            #    if len(seqdict[org]) < 300:
            #        print(org, ':', seqdict[org])

            refseq = seqdict[ref]
            del seqdict[ref]


            if len(seqdict) < 1: # do not align if only one sequence is left
                old = syn
                syn = next(syniter)[1][0]
                continue

            # construct alignment index from the reference
            logger.info("Starting Alignment")
            aligner = mp.Aligner(seq=refseq)#, preset='asm5') 
            alns = {org: align_concatseqs(aligner, seq, chr, mappingtrees[org]) for org, seq in seqdict.items()}
            logger.info(f"None/empty in Alignments: {[org for org in alns if alns[org] is None or alns[org].empty]}")
            print(ref, refseq)
            print(seqdict)


            # run syri
            cwd = '/tmp/' #util.TMPDIR if util.TMPDIR else '/tmp/'
            syris = {org:getsyriout(alns[org], PR='', CWD=cwd) for org in alns if alns[org] is not None and not alns[org].empty}
            # skip regions that were skipped or could not be aligned

            for org in syris:
                print(syris[org].head())
                print(alns[org].head())
                # adjust positions to reference using offsets stored in trees

            # call cross/coresyn, probably won't need to remove overlap
            # incorporate into output
            # has to be sorted, how to do this?
            old = syn
            syn = next(syniter)[1][0]

    except StopIteration as e:
        logger.warning(f'Stopped iteration: {e}')

    return ret

cdef align_concatseqs(aligner, seq, cid, tree):
    """
    Function to align the concatenated sequences as they are and then remap the positions to the positions in the actual genome
    Will split alignments spanning multiple offsets (WIP!).
    """
    m = aligner.map(seq)
    al = deque()
    # traverse alignments
    for h in m:
        al.append([h.r_st+1,
                   h.r_en,
                   h.q_st+1,
                   h.q_en,
                   h.r_en - h.r_st,
                   h.q_en - h.q_st,
                   format((sum([i[0] for i in h.cigar if i[1] == 7]) / (0.01 + sum(
                       [i[0] for i in h.cigar if i[1] in [1, 2, 7, 8]]))) * 100, '.2f'),
                   1,
                   h.strand,
                   h.ctg,
                   cid,
                   "".join(map(lambda x: str(x[0]) + 'MIDNSHP=X'[x[1]], h.cigar))
                   ])

    al = pd.DataFrame(al)
    if al.empty:
        return None
    al[6] = al[6].astype('float')
    al = al.loc[al[6] > 90]
    al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] + al.loc[al[8] == -1, 3]
    al.loc[al[8] == -1, 3] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    al.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    al.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    return al


cdef getsyriout(coords, PR='', CWD='.', N=1, TD=500000, TDOLP=0.8, K=False):
    BRT = 20
    TUC = 1000
    TUP = 0.5
    T = 50

    chrs = list(np.unique(coords.aChr))
    with multiprocessing.Pool(processes=N) as pool:
        pool.map(partial(syri, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, tdgl=TD,tdolp=TDOLP), chrs)

    # Merge output of all chromosomes
    mergeOutputFiles(chrs, CWD, PR)

    #Identify cross-chromosomal events in all chromosomes simultaneously
    getCTX(coords, CWD, chrs, T, BRT, PR, TUC, TUP, N, TD, TDOLP)

    # Recalculate syntenic blocks by considering the blocks introduced by CX events
    outSyn(CWD, T, PR)

    o = getsrtable(CWD, PR)
    if not K:
        for fin in ["synOut.txt", "invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt", "ctxOut.txt", "sv.txt", "notAligned.txt", "snps.txt"]:
            try:
                os.remove(CWD+PR+fin)
            except OSError as e:
                if e.errno != 2:    # 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
                    raise
    return o


