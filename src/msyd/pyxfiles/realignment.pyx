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
from functools import partial

from syri.synsearchFunctions import syri, mergeOutputFiles, outSyn
from syri.tdfunc import getCTX
from syri.writeout import getsrtable

import msyd.util as util
import msyd.classes.cigar as cigar


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
                chrom = syn.ranges_dict[org].chr # chr must always stay the same
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
                        #TODO add Ns in place of known haplotype maybe
                        offset = synrng.end # skip till the end
                        continue

                    #print(l, offset, seq)
                    # add to the intervaltree
                    tree[pos:pos+l] = offset
                    seq += fasta.fetch(region=chrom, start=offset, end=offset+l)
                    offset += l
                    pos += l

                # see if theres any sequence left to realign after processing the crosssyn regions
                l = syn.ranges_dict[org].start - offset
                if l >= MIN_REALIGN_THRESH:
                    tree[pos:pos+l] = offset
                    seq += fasta.fetch(region=chrom, start=offset, end=offset+l)

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
            reftree = mappingtrees[ref]
            del mappingtrees[ref]


            if len(seqdict) < 1: # do not align if only one sequence is left
                old = syn
                syn = next(syniter)[1][0]
                continue

            # construct alignment index from the reference
            logger.info("Starting Alignment")
            aligner = mp.Aligner(seq=refseq, preset='asm5') 
            alns = {org: align_concatseqs(aligner, seq, chrom, reftree, mappingtrees[org]) for org, seq in seqdict.items()}

            # filter out alignments only containing inversions
            for org in alns:
                if alns[org] is not None and all(alns[org].bDir == -1):
                    logger.warning(f"{org} is None in alns or only contains inverted alignments: \n{alns[org]}")
                    alns[org] = None

            logger.info(f"None/empty in Alignments: {[org for org in alns if alns[org] is None]}")
            #print(ref, refseq)
            #print(seqdict)


            # run syri
            cwd = '/tmp/' #util.TMPDIR if util.TMPDIR else '/tmp/'
            syris = {org:getsyriout(alns[org], PR='', CWD=cwd) for org in alns if alns[org] is not None}
            # skip regions that were skipped or could not be aligned

            for org in syris:
                print(syris[org].head())
                print(alns[org].head())

            # call cross/coresyn, probably won't need to remove overlap
            # incorporate into output
            # has to be sorted, how to do this?
            old = syn
            syn = next(syniter)[1][0]

    except StopIteration as e:
        logger.warning(f'Stopped iteration: {e}')

    return ret

cdef align_concatseqs(aligner, seq, cid, reftree, qrytree):
    """
    Function to align the concatenated sequences as they are and then remap the positions to the positions in the actual genome.
    Both sequences should be on the same chromosomes.
    TODO make it split alignments spanning multiple offsets?
    """
    m = aligner.map(seq, extra_flags=0x4000000) # this is the --eqx flag, causing X/= to be added instead of M tags to the CIGAR string
    #print([str(x) for x in m])
    al = deque()
    # traverse alignments
    for h in m:
        rstart: int = h.r_st
        rend: int = h.r_en
        qstart: int = h.q_st
        qend: int = h.q_en
        cg = cigar.cigar_from_bam(h.cigar)
        #print(h.mapq)
        if rstart > rend:
            logger.error(f"Inverted on Reference: {h}")
            continue
            # shouldn't ever occur, TODO maybe handle anyway?

        rstartov = list(reftree[rstart])[0]
        qstartov = list(qrytree[qstart])[0]

        # simply append alignment if there is only one offset
        # as this happens quite often, this should save a lot of time
        if rstartov == list(reftree[rend-1])[0] and qstartov == list(qrytree[qend-1])[0]:
            roff = rstartov.data
            qoff = qstartov.data
            al.append([rstart + roff, rend + roff, qstart + qoff, qend + qoff, rend - rstart, qend - qstart, cg.get_identity()*100, 1 if rstart < rend else -1, h.strand, cid, cid, cg.to_string()])
            continue


        for rint in sorted(reftree[rstart:rend]):
            # subset alignment to this reference offset interval
            qstdel, rcg = cg.get_removed(max(rint.begin - rstart, 0))
            qendel, rcg = rcg.get_removed(max(rend - rint.end, 0), start=False)
            for qint in sorted(qrytree[qstart + qstdel:qend - qendel]):
                # subset to the query offset, respecting the subsetting done so far
                rstdel, qcg = rcg.get_removed(max(qint.begin - qstdel - qstart, 0), ref=False)
                rendel, qcg = qcg.get_removed(max(qend - qint.end - qendel, 0), ref=False, start=False)

                #TODO maybe filter out small alignments here?
                al.append([rint.data + rstdel, rint.data + min(rend, rint.end) - rendel,
                           qint.data + max(qstart, qint.begin), qint.data + min(qend, qint.end),
                           min(rend, rint.end) - rendel - rstdel, min(qend, qint.end) - max(qstart, qint.begin),
                           qcg.get_identity()*100, 1 if rstart < rend else -1, 1 if qstart < qend else -1, cid, cid, qcg.to_string()])

        #for rint in refoffsets:
        #    print("rint:", rint)
        #    # drop from the alignment everything before the current interval
        #    start = max(rint.begin, rstart)
        #    qstartdelta, cg = cg.get_removed(start - rstart)
        #    print("start, rstart, qstartdelta", start, rstart, qstartdelta)
        #    rstart = start

        #    # drop everything after the current interval
        #    end = min(rint.end, rend)
        #    qenddelta, rcg = cg.get_removed(rend - end, start=False)
        #    roffset = rint.data
        #    print("end, rend, qenddelta, roffset", end, rend, qenddelta, roffset)

        #    for qint in sorted(qrytree[qstart + qstartdelta:qend - qenddelta]):
        #        print("qint:", qint)
        #        start = max(qint.begin - qstartdelta, 0)
        #        # drop from the alignment everything before the current interval
        #        rstartdelta, qcg = rcg.get_removed(start, ref=False)
        #        print("start, rstartdelta", start, rstartdelta)

        #        end = min(qint.end, qend)
        #        # drop everything after the current interval
        #        renddelta, qcg = qcg.get_removed(end - qint.end, start=False, ref=False)

        #        # transform coordinates with the offset/alignment information, return 
        #        qoffset = qint.data
        #        rlen = rend - renddelta - rstart - rstartdelta
        #        qlen = qend - qenddelta - qstart - qstartdelta
        #        if rlen < MIN_REALIGN_THRESH or qlen < MIN_REALIGN_THRESH:
        #            print("very small:", rlen, qlen)

        #        al.append([rstart + rstartdelta + roffset, rend - renddelta + roffset,
        #                   qstart + qstartdelta + qoffset, qend - qenddelta + qoffset,
        #                   rlen, qlen, qcg.get_identity()*100, 1, h.strand, cid, cid, qcg.to_string()])



    al = pd.DataFrame(al)
    if al.empty:
        return None
    #print(al[6])
    #al[6] = al[6].astype('float')
    al = al.loc[al[6] > 90]
    al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] + al.loc[al[8] == -1, 3]
    al.loc[al[8] == -1, 3] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    al.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    al.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    #print(al[['aStart', 'aLen', 'bStart', 'bLen', 'iden']])

    #TODO use tree to remap!

    return None if al.empty else al

cdef subset_ref_offset(rstart, rend, qstart, qend, cg, interval):
    """DEPRECATED
    Takes an alignment and an interval from the intervaltree, returns the part of the alignment that is in the interval on the reference with the offset incorporated
    """
    start = max(interval.start, rstart)
    # drop from the alignment everything before the current interval
    qstartdelta, curcg = cg.get_removed(start - rstart)

    # drop everything after the current interval
    end = min(interval.end, rend)
    qenddelta, retcg = curcg.get_removed(rend - end, start=False)

    # transform coordinates with the offset/alignment information, return 
    offset = interval.data
    return (start + offset, end + offset, qstart + qstartdelta, qend - qenddelta, retcg)

cdef subset_qry_offset(rstart, rend, qstart, qend, cg, interval):
    """DEPRECATED
    Takes an alignment and an interval from the intervaltree, returns the part of the alignment that is in the interval on the query with the offset incorporated
    """
    start = max(interval.start, qstart)
    # drop from the alignment everything before the current interval
    rstartdelta, curcg = cg.get_removed(start - qstart, ref=False)

    end = min(interval.end, qend)
    # drop everything after the current interval
    renddelta, retcg = curcg.get_removed(end - interval.end, start=False, ref=False)

    # transform coordinates with the offset/alignment information, return 
    offset = interval.data
    return (rstart + rstartdelta, rend + renddelta, start + offset, end + offset, retcg)


cdef getsyriout(coords, PR='', CWD='.', N=1, TD=500000, TDOLP=0.8, K=False):
    BRT = 20
    TUC = 1000
    TUP = 0.5
    T = 50
    invgl = 1000000

    #chrs = list(np.unique(coords.aChr))
    assert(len(list(np.unique(coords.aChr))) == 1)
    print(coords[['aStart', 'aEnd', 'aLen', 'bStart', 'bEnd', 'bLen', 'iden', 'aDir', 'bDir']])#, 'cigar']])
    chrom = list(coords.aChr)[0] # there should only ever be one chr anyway
    syri(chrom, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD, tdolp=TDOLP)

    #with multiprocessing.Pool(processes=N) as pool:
    #    pool.map(partial(syri, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD,tdolp=TDOLP), chrs)


    #TODO if runtime a problem: redo syri call to only call synteny => maybe configurable?
    # Merge output of all chromosomes – still necessary for some reason
    mergeOutputFiles([chrom], CWD, PR)

    #Identify cross-chromosomal events in all chromosomes simultaneously
    getCTX(coords, CWD, [chrom], T, BRT, PR, TUC, TUP, N, TD, TDOLP)

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


