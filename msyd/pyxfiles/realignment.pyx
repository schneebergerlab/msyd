#%%cython
#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import sys

import pandas as pd
import numpy as np
import mappy as mp
import pysam
import intervaltree
from datetime import datetime
from multiprocessing import Pool
import logging
from collections import deque, defaultdict
import os
from functools import partial
from io import StringIO

cimport libc.stdio as cio
cimport posix.unistd as unistd

# I added these lines to hide all of the INFO logs from syri. If those are required then these lines can be removed
logging.getLogger('syri').setLevel(logging.WARNING)
logging.getLogger('getCTX').setLevel(logging.WARNING)
from syri.synsearchFunctions import syri, mergeOutputFiles, outSyn, apply_TS, alignmentBlock, getSynPath
from syri.tdfunc import getCTX
from syri.writeout import getsrtable

import msyd.scripts.util as util
import msyd.cigar as cigar
import msyd.pansyn as pansyn
import msyd.io as io


cdef int _MIN_REALIGN_THRESH = 100 # min length to realign regions
cdef int _MAX_REALIGN = 0 # max number of haplotypes to realign to
cdef int _NULL_CNT = 20 # number of separators to use between blocks during alignment

logger = util.CustomFormatter.getlogger(__name__)

# <editor-fold desc='Support functions for realign'>
cpdef construct_mappingtrees(crosssyns, old, syn):
# def construct_mappingtrees(crosssyns, old, syn):
    """
    Makes a dictionary containing an intervaltree with an offset mapping for each org containing enough non-aligned sequence to realign.
    Crosssyns need to be sorted by position on reference.
    For each tree, the sequence in genome `org` at position `tree[pos].data - tree[pos].begin + pos` corresponds to the position `pos` in the synthetic query sequence.
    """
    mappingtrees = defaultdict(intervaltree.IntervalTree)
    posdict = defaultdict(int) # stores the current position in each org
    offsetdict = {org:rng.end for org, rng in old.ranges_dict.items()} # stores the current offset in each org

    for reforg in crosssyns:
        for crosssyn in crosssyns[reforg]:
            # iterate through all pansyns found so far by realignment
            for org, rng in crosssyn.ranges_dict.items():
                #print(f"{offsetdict[org]}, {posdict[org]}, {rng}, {mappingtrees[org]}")
                l = rng.start - offsetdict[org] # len of the region to be added
                if l < 0:
                    # print('ERROR')# improper sorting – skip
                    continue

                # check if this interval would be redundant
                prev = list(mappingtrees[org][posdict[org]-1]) # will be empty if tree is empty
                if prev and posdict[org] + prev[0].data == offsetdict[org]:
                    # extend the previous interval instead
                    del mappingtrees[org][posdict[org]-1]
                    posdict[org] += l
                    mappingtrees[org][prev[0].begin:posdict[org]] = prev[0].data

                elif l > _MIN_REALIGN_THRESH: # otherwise add to the tree if it's large enough
                    mappingtrees[org][posdict[org]:posdict[org]+l] = offsetdict[org]
                    posdict[org] += l + _NULL_CNT # add the spacer length to the start of the next index

                # all up to the end of this region has been added
                offsetdict[org] = rng.end

    # see if there's any sequence left to realign after processing the crosssyn regions
    for org, offset in offsetdict.items():
        l = syn.ranges_dict[org].start - offset
        if l >= _MIN_REALIGN_THRESH:
            mappingtrees[org][posdict[org]:posdict[org]+l] = offset
    return mappingtrees
# END

cpdef align_concatseqs(seq, qcid, qrytree, refseq, preset, rcid, reftree, aligner=None):
# def align_concatseqs(seq, qcid, qrytree, refseq, preset, rcid, reftree):
    """
    Function to align the concatenated sequences as they are and then remap the positions to the positions in the actual genome.
    Both sequences should be on the same chromosomes.
    Splits alignments that span multiple offsets into one alignment per offset
    """
    # Parse aligner from parent function when not using multiprocessing.Pool. When using Pool, define aligner here
    if aligner is None:
        aligner = mp.Aligner(seq=refseq, preset=preset)
    # logger.debug('Start align_concatseqs')
    m = aligner.map(seq, extra_flags=0x4000000) # this is the --eqx flag, causing X/= to be added instead of M tags to the CIGAR string
    logger.debug(f'Minimap2 alignment done.')
    #print([str(x) for x in m])
    al = deque()
    # traverse alignments
    logger.debug('Traversing alignments')
    for h in m:
        # print('a')
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
        # print(f'reftree: {reftree}, qrytree: {qrytree}, rend: {rend}, qend: {qend}')

        if rstartov == list(reftree[rend-1])[0] and qstartov == list(qrytree[qend-1])[0]:
            roff = rstartov.data
            qoff = qstartov.data
            al.append([rstart + roff, rend + roff, qstart + qoff, qend + qoff, rend - rstart, qend - qstart, cg.get_identity()*100, 1 if rstart < rend else -1, h.strand, rcid, qcid, cg.to_string()])
            continue
        # print('b')

        for rint in sorted(reftree[rstart:rend]):
            # subset alignment to this reference offset interval
            qstdel, rcg = cg.get_removed(max(rint.begin - rstart, 0))
            qendel, rcg = rcg.get_removed(max(rend - rint.end, 0), start=False)
            for qint in sorted(qrytree[qstart + qstdel:qend - qendel]):
                # subset to the query offset, respecting the subsetting done so far
                rstdel, qcg = rcg.get_removed(max(qint.begin - qstdel - qstart, 0), ref=False)
                rendel, qcg = qcg.get_removed(max(qend - qint.end - qendel, 0), ref=False, start=False)

                #TODO maybe filter out small alignments here?
                #print("r:", rint.data, rstart, rend, rint.begin, rint.end, rendel, rstdel, qcg.get_len(ref=True))
                #print("q:", qint.data, qstart, qend, qint.begin, qint.end, qendel, qstdel, qcg.get_len(ref=False))
                al.append([rint.data + rstdel, rint.data + min(rend, rint.end) - rendel - max(rint.begin - rstart, 0),
                           qint.data + max(qstart, qint.begin), qint.data + min(qend, qint.end),
                           min(rend, rint.end) - rendel - rstdel - max(rint.begin - rstart, 0), min(qend, qint.end) - max(qstart, qint.begin),
                           qcg.get_identity()*100, 1 if rstart < rend else -1, 1 if qstart < qend else -1, rcid, qcid, qcg.to_string()])
        # print('c')
    logger.debug('Alignments traversed')
    # print('d')
    al = pd.DataFrame(al)
    if al.empty:
        return None
    #print(al[6])
    #al[6] = al[6].astype('float')
    al = al.loc[al[6] > 90] # TODO: Alignment identity filter. This filter is not mandatory and the user might opt to remove this
    al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] + al.loc[al[8] == -1, 3]
    al.loc[al[8] == -1, 3] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    al.loc[al[8] == -1, 2] = al.loc[al[8] == -1, 2] - al.loc[al[8] == -1, 3]
    al.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    al.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    #print(al[['aStart', 'aLen', 'bStart', 'bLen', 'iden']])
    # print('e')
    #TODO use tree to remap!
    return None if al.empty else al

# </editor-fold>

cpdef realign(df, qrynames, fastas, MIN_REALIGN_THRESH=None, MAX_REALIGN=None, NULL_CNT=None, mp_preset='asm5', ncores=1):
    if MIN_REALIGN_THRESH is not None and MIN_REALIGN_THRESH >= 0:
        global _MIN_REALIGN_THRESH
        _MIN_REALIGN_THRESH = int(MIN_REALIGN_THRESH)
    if MAX_REALIGN is not None and MAX_REALIGN >= 0:
        global _MAX_REALIGN
        _MAX_REALIGN = int(MAX_REALIGN)
    if NULL_CNT is not None and NULL_CNT >= 0:
        global _NULL_CNT
        _NULL_CNT = int(NULL_CNT)

    # return process_gaps(df, qrynames, fastas, mp_preset=mp_preset, ncores=ncores, cwd=cwd)

# cpdef process_gaps(df, qrynames, fastas, mp_preset, ncores, cwd):
    """
    Function to find gaps between two coresyn regions and realign them to a new reference.
    Discovers all crosssynteny.
    
    :arguments: A DataFrame with core and crosssyn regions called by find_pansyn and the sample genomes and names. `mp_preset` designates which minimap2 alignment preset to use.
    :returns: A DataFrame with the added non-reference crosssynteny
    """
    # init stuff
    ret = deque()#pd.DataFrame()
    n = len(qrynames)
    if not n == len(fastas):
        logger.error(f"More/less query names than fastas passed to process_gaps: {qrynames}, {fastas}")
        raise ValueError("Wrong number of fastas!")

    # create single-byte unique identifiers for each query that aren't ACTGX
    # this is required to fill the alignments, as minimap2 will align non-bases just like bases
    # this chooses the first character in each sample name that is both free and legal to use.
    #TODO this is a pretty broken system, and fails if using more than 21 samples (assuming english input)
    # either find a way to make minimap work or extend to use multiple characters if required (super inefficient though)
    filler_dict = {org: '' for org in qrynames}
    if _NULL_CNT > 0:
        forbidden = set(['A', 'C', 'G', 'T', 'X'])
        for org in qrynames:
            for ch in org:
                ch = ch.upper()
                if ch not in forbidden and not ch in filler_dict.values():
                    filler_dict[org] = ch
                    break

            if not filler_dict[org]:
                logger.error(f"Unable to find unique characters for every org in {qrynames}. Final mapping {filler_dict}. This could be because there are more than 21 samples. Try calling with NULL_CNT set to 0.")
                raise ValueError("Unable to find unique characters for every org!")

    logger.info(filler_dict)


    # load fasta files
    fafin = {qrynames[i]: pysam.FastaFile(fastas[i]) for i in range(len(qrynames))}

    # iterate through each gap between coresyn blocks
    # call the alignment/ functionality and merge   

    syniter = df.iterrows()
    try:
        syn = next(syniter)[1][0]
        # skip to first core
        while syn.get_degree() < n:
            ret.append(syn)
            syn = next(syniter)[1][0]
            # print(vars(syn))
        old = syn
        # TODO: Misses the crosssyn before the first coresyn region? This would become if there are no or very few coresyn regions
        # leon: yes, the crosssyn before the first and after the last coresyn of a chromosome will be missed.
        # my thinking was that this would correspond to the highly polymorphic tails of the chromosomes
        # this seems to work out in Ampril (the first coresyn starts at <2 kb in my local test dataset),
        # but might cause problems in datasets with very little coresyn
        crosssyns = dict()
        # CNT = 0
        while True:
            # CNT += 1
            # print(CNT)
            # if CNT == 100:
            #     break
            # find block between two coresyn regions
            crosssyns = dict()

            # store crosssyn-regions, to be added once we know if we need to realign
            refcrosssyns = []
            while syn.get_degree() < n:
                refcrosssyns.append(syn)
                syn = next(syniter)[1][0]

            crosssyns[old.ref.org] = refcrosssyns

            # syn must be core now

            # start and end of the non-ref region, on the reference
            end = syn.ref.start
            start = old.ref.end

            # if there is not enough novel sequence on any organism to realign, skip this realignment preemptively
            if end - start < _MIN_REALIGN_THRESH:
                if all([syn.ranges_dict[org].start - old.ranges_dict[org].end < _MIN_REALIGN_THRESH \
                        for org in syn.ranges_dict]):
                    ret.append(syn)
                    old = syn
                    syn = next(syniter)[1][0]
                    continue

            # between chromosomes, there isn't a single gap
            #print(f"Chrs: {syn.ref.chr}, {old.ref.chr}")
            if syn.ref.chr != old.ref.chr:
                logger.error("Chr case found: {syn.ref}, {old.ref}")
                old = syn
                ret.append(syn)
                syn = next(syniter)[1][0]
                continue

            # Block has been extracted and is long enough;
            # extract appropriate sequences, respecting crossyn

            #print(old, syn, syn.ref.start - old.ref.end)
            #print(crosssyns)

            ##TODO parallelise everything after this, enable nogil


            # construct a mapping tree and concatenate all the sequences contained within
            mappingtrees = construct_mappingtrees(crosssyns, old, syn)
            seqdict = {org: (filler_dict[org]*_NULL_CNT).join([
                fafin[org].fetch(region = syn.ranges_dict[org].chr,
                                 start = interval.data - (ind*_NULL_CNT), # subtract the spacers before this point
                                 end = interval.data + interval.end - interval.begin - (ind*_NULL_CNT))
                for ind, interval in enumerate(sorted(mappingtrees[org]))]) # in theory intervaltrees should sort itself, but just in case
                for org in mappingtrees}

            #if any([len(mappingtrees[org]) > 3 for org in mappingtrees]):
            #    logger.info(mappingtrees)
            #    logger.info(seqdict)

            if not seqdict: # if all sequences have been discarded, skip realignment
                #logger.info("Not aligning, not enough non-reference sequence found!")
                old = syn
                ret.append(syn)
                syn = next(syniter)[1][0]
                continue

            # TODO: Parse file containing the centromere coordinate. Check if the selected region is centromeric, skip re-alignment if it is.

            while len(seqdict) > 2: # realign until there is only one sequence left

                # stop realignment if we have already found _MAX_REALIGN haplotypes
                if _MAX_REALIGN > 0 and len(crosssyns) > _MAX_REALIGN:
                    break

                # choose a reference as the sample containing the most non-crosssynteny
                # ref = max(map(lambda x: (len(x[1]), x[0]), seqdict.items()))[1]
                ref = max([(len(v), k) for k,v in seqdict.items()])[1]

                #print('ref:', ref)
                #print('On ref:', syn.ref.chr, start, end, end - start)
                #print({org:len(seq) for org, seq in seqdict.items()})
                #        print(org, ':', seqdict[org])

                refseq = seqdict[ref]
                del seqdict[ref]
                reftree = mappingtrees[ref]
                del mappingtrees[ref]

                # construct alignment index from the reference
                logger.debug(f"Starting Alignment. Left core: {old.ref}. Right core: {syn.ref}")
                # aligner = mp.Aligner(seq=refseq, preset=mp_preset)
                # print('start alignment', datetime.now())
                # alns = {}
                # TODO: The alignment step is a major performance bottleneck, specially when aligning centromeric regions. If the expected memory load is not high, then we can easily parallelise align_concatseqs using multiprocessing.Pool. Here, I have implemented it hoping that it should not be a problem. If at some point, we observe that the memory footprint increases significantly, then we might need to revert it back.
                # print(1)
                if syn.ref.start - old.ref.end > 100000:
                    # print('par')
                    alignargs = [[seqdict[org], syn.ranges_dict[org].chr, mappingtrees[org]] for org in seqdict.keys()]
                    with Pool(processes=ncores) as pool:
                        # pool.starmap(partial(foo, d='x'), alignargs)
                        alns = pool.starmap(partial(align_concatseqs, refseq=refseq, preset=mp_preset, rcid=syn.ref.chr, reftree=reftree, aligner=None), alignargs)
                    alns = dict(zip(list(seqdict.keys()), alns))
                else:
                    aligner = mp.Aligner(seq=refseq, preset=mp_preset)
                    alns = dict()
                    # print('seq')
                    for org, seq in seqdict.items():
                        # print(org, datetime.now())
                        # TODO: Currently (12.03.2024), this seems to be the most time-consuming step
                        alns[org] = align_concatseqs(seq, syn.ranges_dict[org].chr, mappingtrees[org], refseq, mp_preset, syn.ref.chr, reftree, aligner=aligner)


                # filter out alignments only containing inversions
                for org in alns:
                    if alns[org] is not None and all(alns[org].bDir == -1):
                        logger.warning(f"{org} in alns only contains inverted alignments: \n{alns[org]}")
                        alns[org] = None

                #logger.info(f"None/empty in Alignments: {[org for org in alns if alns[org] is None]}")
                #print(ref, refseq)
                #print(seqdict)

                # run syri
                logger.debug("Running syri")

                # TODO MG: Replaced getsyriout with the synteny identification method from syri. Consider parallelizing this because for repetitive regions this loop would be expensive as well.
                syris = syri_get_syntenic(alns)

                for org in syris:
                    if syris[org] is not None:
                        #print("===", org, mappingtrees[org][0], "against", ref, reftree[0], reftree[-1], "===")
                        #print(syris[org])

                        alns[org].columns = ["astart", "aend", "bstart", "bend", "alen", "blen", "iden", "adir", "bdir", "achr", "bchr", 'cg']
                        #print(alns[org][['astart', 'aend', 'alen', 'bstart', 'bend', 'blen', 'bdir', 'iden']])
                        #print(mappingtrees[org])
                syns = [pansyn.match_synal(
                            io.extract_syri_regions(syris[org], reforg=ref, qryorg=org, anns=["SYNAL"]),
                            alns[org])#, ref=ref)
                        for org in syris if syris[org] is not None]


                if len(syns) == 0:
                    continue

                # syns should be sorted
                # TODO: MG Question: what is the use of this function?
                # leon: I've found that sometimes, syri calls would overlap a bit,
                # which the pansyn identification doesn't like
                # this checks this isn't the case and corrects it in case it finds it
                # should probably not be strictly necessary here as we call syri ourselves
                # but doesn't hurt to check either

                # print(syns)
                pansyns = pansyn.reduce_find_overlaps(syns, cores=1)
                # print(6)
                # no need to recalculate the tree if no pansynteny was found
                if pansyns is None or pansyns.empty:
                    continue

                # Add all crosssyns with alphabetical sorting by reference name
                crosssyns[ref] = [psyn[1][0] for psyn in pansyns.iterrows()]
                added = sum([len(x.ref) for x in crosssyns[ref]])

                logger.info(f"Realigned {old.ref.chr}:{old.ref.end}-{syn.ref.start} (len {util.siprefix(syn.ref.start - old.ref.end)}) to {ref}. Found {util.siprefix(added)} (avg {util.siprefix(added/len(crosssyns))}) of cross-synteny.")

                ## recalculate mappingtrees from current crosssyns to remove newly found cross synteny
                # TODO maybe in future directly remove, might be more efficient
                mappingtrees = construct_mappingtrees(crosssyns, old, syn)
                # remove all orgs that have already been used as a reference
                for reforg in crosssyns:
                    if reforg in mappingtrees:
                        del mappingtrees[reforg]

                seqdict = {org:(filler_dict[org]*_NULL_CNT).join([
                    fafin[org].fetch(region = syn.ranges_dict[org].chr,
                                     start = interval.data - (ind*_NULL_CNT), # subtract the spacers before this point
                                     end = interval.data + interval.end - interval.begin - (ind*_NULL_CNT))
                    for ind, interval in enumerate(sorted(mappingtrees[org]))])
                    for org in mappingtrees}
                if not seqdict: # if all sequences have been discarded, finish realignment
                    break

            # incorporate into output DF, sorted alphabetically by ref name
            # does nothing if no crossyn was found
            for org in sorted(crosssyns.keys()):
                ret.extend(crosssyns[org])

            # continue checking the next coresyn gap
            old = syn
            ret.append(syn)
            syn = next(syniter)[1][0]
            #print(old, syn)
    except StopIteration as e:
        logger.warning(f'Stopped iteration: {e}')
    return pd.DataFrame(list(ret))
# END


cdef syri_get_syntenic(alns):
    # Synteny call parameters
    BRT = 20
    TUC = 1000
    TUP = 0.5
    T = 50
    invgl = 1000000

    syris = {}

    #NOTE: for large regions, it might make sense to parallelize the syri call
    # per organism, turning the for loop below into a parallelized map
    # probably only worth it with no_gil, though
    for org in alns:
        if alns[org] is None: continue
        coords = alns[org]
        try:
            assert coords.aChr.nunique() == 1
            assert coords.bChr.nunique() == 1
        except AssertionError:
            logger.error(
                f"Incorrect coords. More than one chromosome parsed. Ref chromosomes: {coords.aChr}. Qry chromosomes: {coords.bChr}")
        # NOTE: syri requires that the coords table have same chromosome IDs for homologous chromosomes. When, the coords have different chromosome IDs, then manipulate the chroms IDs here
        chromr = list(coords.aChr)[0]  # there should only ever be one chr anyway
        chromq = list(coords.bChr)[0]
        samechrids = chromr == chromq
        if not samechrids:
            coords.bChr.replace(chromq, chromr, inplace=True)
        chromo = chromr
        coordsData = coords[(coords.aChr == chromo) & (coords.bChr == chromo) & (coords.bDir == 1)]
        syndf = apply_TS(coordsData.aStart.values, coordsData.aEnd.values, coordsData.bStart.values,
                      coordsData.bEnd.values, T)

        blocks = [alignmentBlock(i, syndf[i], coordsData.iloc[i]) for i in syndf.keys()]
        for block in blocks:
            i = 0
            while i < len(block.children):
                block.children = list(set(block.children) - set(blocks[block.children[i]].children))
                i += 1
            block.children.sort()
            for child in block.children:
                blocks[child].addParent(block.id)
            scores = [blocks[parent].score for parent in block.parents]
            if len(scores) > 0:
                block.bestParent(block.parents[scores.index(max(scores))], max(scores))

        synPath = getSynPath(blocks)
        synData = coordsData.iloc[synPath].copy()

        if not samechrids:
            synData.bChr.replace(chromr, chromq, inplace=True)

        synData.columns = list(map(str.lower, synData.columns))
        synData[['aseq', 'bseq', 'id', 'parent', 'dupclass']] = '-'     # Setting `id` and `parent` as `-`. This does not fit normal syri output but here it should inconsequential as these columns are (probably) not used anyway

        synData['vartype'] = 'SYNAL'
        synData = synData[['achr', 'astart', 'aend', 'aseq', 'bseq', 'bchr', 'bstart', 'bend', 'id', 'parent', 'vartype', 'dupclass']]
        syris[org] = synData
        # print(synData.columns)
    # skip regions that were skipped or could not be aligned, or only contain inverted alignments

    return syris



################################################# DEPRECATED ###########################################################
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




# TODO: Make parameters adjustable. Also, now (12.03.24) this could be deprecated.
cpdef getsyriout(coords, PR='', CWD='.', N=1, TD=500000, TDOLP=0.8, K=False, redir_stderr=False):
    """DEPRECATED"""
    BRT = 20
    TUC = 1000
    TUP = 0.5
    T = 50
    invgl = 1000000

    #assert(len(list(np.unique(coords.aChr))) == 1)
    try:
        assert coords.aChr.nunique() == 1
        assert coords.bChr.nunique() == 1
    except AssertionError:
        logger.error(f"Incorrect coords. More than one chromosome parsed. Ref chromosomes: {coords.aChr}. Qry chromosomes: {coords.bChr}")

    cdef int oldstderr = -1
    if redir_stderr:
        #cio.fclose(cio.stderr)
        #cio.stderr = cio.freopen(bytes(f"{CWD}/stderr", encoding='utf8'), "w", cio.stderr)
        oldstderr = unistd.dup(unistd.STDERR_FILENO)
        cio.freopen(bytes(f"{CWD}/stderr", encoding='utf8'), "w", cio.stderr)

    # NOTE: syri requires that the coords table have same chromosome IDs for homologous chromosomes. When, the coords have different chromosome IDs, then manipulate the chroms IDs here
    chromr = list(coords.aChr)[0] # there should only ever be one chr anyway
    chromq = list(coords.bChr)[0]
    samechrids = chromr == chromq
    if not samechrids:
        coords.bChr.replace(chromq, chromr, inplace=True)
    chrom = chromr
    # handle errors by return value; allows only showing output if there is a problem
    # python errors coming after an error here will have normal stderr
    # try:
    #     # TODO: this function expects that the reference and query chromsome would have the same id. If that is not the case (pre-processing not done),then this function always crashes
    #     syri(chrom, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD,tdolp=TDOLP)
    # except ValueError:
    #     print(coords[['aStart', 'aEnd', 'aLen', 'bStart', 'bEnd', 'bLen', 'iden', 'aDir', 'bDir']])
    #     return None

    if syri(chrom, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD, tdolp=TDOLP) == -1:
        if redir_stderr:
            logger.error("Redirecting stderr to console again")
            #cio.fclose(cio.stderr)
            #cio.stderr = oldstderr
            unistd.close(unistd.STDERR_FILENO)
            unistd.dup2(oldstderr, unistd.STDERR_FILENO)
        logger.error("syri call failed on input:")
        print(coords[['aStart', 'aEnd', 'aLen', 'bStart', 'bEnd', 'bLen', 'iden', 'aDir', 'bDir']])
        if redir_stderr:
            logger.error(f"syri stderr in '{CWD}/stderr'")
        return None

    #with multiprocessing.Pool(processes=N) as pool:
    #    pool.map(partial(syri, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD,tdolp=TDOLP), chrs)

    #TODO if runtime a problem: redo syri call to only call synteny => maybe configurable?
    # Merge output of all chromosomes – still necessary for some reason
    mergeOutputFiles([chrom], CWD, PR)

    #TODO: Maybe not requires and can be removed?
    # leon: outSyn fails if this isn't called, that's why it's left in there
    # but yes, this step should be unnecessary
    # In general, I think the syri calling should be done more elegantly --
    # writing to and then reading from files is quite inefficient, especially
    # for short realignments

    #Identify cross-chromosomal events in all chromosomes simultaneously
    getCTX(coords, CWD, [chrom], T, BRT, PR, TUC, TUP, N, TD, TDOLP)

    # Recalculate syntenic blocks by considering the blocks introduced by CX events
    outSyn(CWD, T, PR)
    o = getsrtable(CWD, PR)
    if not samechrids:
        o.bchr.replace(chromr, chromq, inplace=True)

    if redir_stderr:
        #cio.fclose(cio.stderr)
        #cio.stderr = oldstderr
        unistd.close(unistd.STDERR_FILENO)
        unistd.dup2(oldstderr, unistd.STDERR_FILENO)

    if not K:
        for fin in ["synOut.txt", "invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt", "ctxOut.txt", "sv.txt", "notAligned.txt", "snps.txt"]:
            try:
                os.remove(CWD+PR+fin)
            except OSError as e:
                if e.errno != 2:    # 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
                    raise
    return o
# END

