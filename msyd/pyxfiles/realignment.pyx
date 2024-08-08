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
import logging
import os
from collections import deque, defaultdict
from functools import partial
from multiprocessing import Pool

from intervaltree import IntervalTree, Interval

# I added these lines to hide all of the INFO logs from syri. If those are required then these lines can be removed
logging.getLogger('syri').setLevel(logging.WARNING)
logging.getLogger('getCTX').setLevel(logging.WARNING)

from syri.synsearchFunctions import syri, mergeOutputFiles, outSyn, apply_TS, alignmentBlock, getSynPath
from syri.tdfunc import getCTX
from syri.writeout import getsrtable

import msyd.util as util
import msyd.cigar as cigar
import msyd.pansyn as pansyn
import msyd.io as io


cdef int _MIN_REALIGN_THRESH = 100 # min length to realign regions
cdef int _MAX_REALIGN = 0 # max number of haplotypes to realign to; set to 0 to realign without limit
cdef int _NULL_CNT = 100 # number of separators to use between blocks during alignment

logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)

cpdef listdict_to_mts(lists):
    """
    Transforms a list of alignments for each organism into a dictionary of intervaltrees mapping an offset in the alignment reference to each position in the sample genome.
    """
    ret = dict()
    posdict = defaultdict(int) # stores the current position in each org
    for org, lst in lists.items():
        tree = IntervalTree()
        for offset, length in lst:
            if length > _MIN_REALIGN_THRESH: # filter again for sufficient len
                tree[posdict[org]:posdict[org] + length] = offset
                posdict[org] += length + _NULL_CNT # add interval + spacer

        #if len(tree) > 0:
        ret[org] = tree
            
    #return {org: IntervalTree(functools.reduce(Interval())
    #        for org, lst in lists.items()}
    return ret

# <editor-fold desc='Support functions for realign'>
cpdef construct_mts(merasyns, old, syn):
# def construct_mappingtrees(merasyns, old, syn):
    """
    Makes a dictionary containing an intervaltree with an offset mapping for each org containing enough non-aligned sequence to realign.
    Crosssyns need to be sorted by position on reference.
    For each tree, the sequence in genome `org` at position `tree[pos].data - tree[pos].begin + pos` corresponds to the position `pos` in the synthetic query sequence.
    """
    listdict = defaultdict(list)
    offsetdict = {org:rng.end for org, rng in old.ranges_dict.items()} # stores the current offset in each org

    for merasyn in merasyns:
        # iterate through all pansyns found so far
        for org, rng in merasyn.ranges_dict.items():
            #print(f"{offsetdict[org]}, {posdict[org]}, {rng}, {mappingtrees[org]}")
            l = rng.start - offsetdict[org] # len of the region to be added
            if l < 0:
                logger.error(f"improper sorting: {rng.start} < {offsetdict[org]}") # improper sorting – skip
                continue

            # check if this interval would be redundant
            if org in listdict:
                prev = listdict[org][-1]
                #print(prev[0] + prev[1], offsetdict[org])
                if prev[0] + prev[1] == offsetdict[org]: 
                    # check if offset + len matches the current offset; extend prev interval instead
                    # no +1 because the end is not inclusive
                    prev[1] += l

            if l > _MIN_REALIGN_THRESH: # otherwise add to the tree if it's large enough
                listdict[org].append( (offsetdict[org], l) )

            # all up to the end of this region has been added
            offsetdict[org] = rng.end

    # see if there's any sequence left to realign after processing the merasyn regions
    for org, offset in offsetdict.items():
        l = syn.ranges_dict[org].start - offset
        if l >= _MIN_REALIGN_THRESH:
            listdict[org].append( (offset, l) )

    for org in old.ranges_dict:
        if old.ranges_dict[org].end > syn.ranges_dict[org].start:
            logger.error(f"{org}: End ({old.ranges_dict[org].end}) after start ({syn.ranges_dict[org].start})! {old} (check: {old.check()}), {syn} (check: {syn.check()}).")
            #logger.debug(f"CIGARS of above error: {old.cigars_dict[org].to_string()}, {syn.cigars_dict[org].to_string()}")

    return listdict_to_mts(listdict)
# END

cpdef subtract_mts(mappingtrees, merasyns):
    """
    Takes a dict containing an `Intervaltree` with offsets for each organism to be realigned (as produced by `construct_mts`), and returns new mappingtrees with regions covered by pansyn objects in `merasyns` subtracted.
    Used to remove merasynteny found during realignment from the mappingtrees, to do further realignment.
    """
    # core iteration is the same as in construct_mts, except we don't need to store the offsets
    curdict = {org:list(mappingtrees[org][0])[0] for org in mappingtrees if len(mappingtrees[org]) > 0} # stores the current interval in each org
    listdict = defaultdict(list) # used to construct the output mappingtrees
    # these need to be reconstructed to take care of handling the separator intervals

    for merasyn in merasyns:
        for org, rng in merasyn.ranges_dict.items():
            if not org in curdict or curdict[org] is None: # skip if there are no more intervals to process for this org
                continue
            curint = curdict[org]
            orglist = listdict[org]
            #TODO idea: set curdict[org] to None to indicate termination?

            # skip to first interval overlapping this merasyn
            # there should be one of thesejfor any alignment
            while curint.end - curint.begin + curint.data < rng.start:
                orglist.append( (curint.data, curint.end - curint.begin) )
                if curint.end + _NULL_CNT + 1 in mappingtrees[org]:
                    curint = list(mappingtrees[org][curint.end + _NULL_CNT +1 ])[0]
                else: # there is no offset after this
                    break
            else: # to skip to next org in big for loop
                #logger.warning("Skipped out of int iteration loop!")
                curdict[org] = None
                continue

            # there shouldn't ever be an alignment spanning beyond one offset, as genuinely adjacent offsets are compressed
            # throw an error if this isn't the case
            #TODO maybe handle improperly compressed offsets? perhaps by ratcheting over the end as welll as in the while loop above
            #assert rng.end <= curint.end - curint.begin + curint.data, "Synteny in a spacer offset detected! An alignment went into the separator. Most likely, something went wrong during alignment."
            if rng.end > curint.end - curint.begin + curint.data:
                logger.debug(f"{rng.end}, {curint.end - curint.begin + curint.data}")
                logger.warning("Synteny in a spacer detected! An alignment went into the separator. Most likely, something went wrong during the alignment call.")
            

            # there was no interval overlapping this merasyn anyway, we don't need to subtract anything
            if curint.data > rng.end:
                curdict[org] = curint
                continue
            
            ## from here, rng is fully within curint

            # remove overlap from start, add as separate offset if large enough
            l = rng.start - curint.data
            if l > _MIN_REALIGN_THRESH:
                orglist.append( (curint.data, l) )

            # set curint to what remains after removing this merasyn if large enough
            l = curint.data + curint.end - curint.begin - rng.end
            if l > _MIN_REALIGN_THRESH:
                curdict[org] = Interval(curint.end - l, curint.end, rng.end)
                #orglist.append( (rng.end, l) )
            else:
                # otherwise skip to next interval
                if curint.end + _NULL_CNT + 1 in mappingtrees[org]:
                    curdict[org] = list(mappingtrees[org][curint.end + _NULL_CNT + 1])[0]
                else:
                    curdict[org] = None
                    continue

    return listdict_to_mts(listdict)


cdef get_aligner(seq, preset, ns=True):
    #aligner = mp.Aligner(seq=refseq, preset=preset)

    # set --score-N parameter to 10
    #aligner.map_opt.sc_ambi = 10

    #aligner = mp.Aligner(seq=refseq, preset=preset, scoring=[1, 19, 39, 81, 3, 1, 10])
    # using values from https://github.com/lh3/minimap2/blob/9b0ff2418c298f5b5d0df12b137896e5c3fb0ef4/options.c#L134
    # https://github.com/lh3/minimap2/issues/155
    # https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/python/README.rst?plain=1#L93
    # values from the manpage, under presets -> asm5
    #-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50
    
    #aligner = mp.Aligner(seq=seq, preset=preset, scoring=[1, 19, 39, 81, 39, 81, 100]) if ns else mp.Aligner(seq=seq, preset=preset)

    aligner = mp.Aligner(seq=seq, preset=preset, sc_ambi=15, max_chain_skip=255) # requires a patched version of minimap2; TODO make PR to get that merged


    return aligner

cpdef align_concatseqs(seq, qcid, qrytree, refseq, preset, rcid, reftree, aligner=None):
# def align_concatseqs(seq, qcid, qrytree, refseq, preset, rcid, reftree):
    """
    Function to align the concatenated sequences as they are and then remap the positions to the positions in the actual genome.
    Both sequences should be on the same chromosomes.
    Splits alignments that span multiple offsets into one alignment per offset
    """
    # Parse aligner from parent function when not using multiprocessing.Pool. When using Pool, define aligner here
    if aligner is None:
        aligner = get_aligner(refseq, preset)

    # logger.debug('Start align_concatseqs')
    m = aligner.map(seq, extra_flags=0x4000000) # this is the --eqx flag, causing X/= to be added instead of M tags to the CIGAR string
    #logger.debug(f'Minimap2 alignment done.')
    #print([str(x) for x in m])
    alns = deque()
    # traverse alignments
    #logger.debug('Traversing alignments')
    for h in m:
        # print('a')
        rstart: int = h.r_st
        rend: int = h.r_en -1 # transform to inclusive range
        qstart: int = h.q_st
        qend: int = h.q_en -1 # transform to inclusive range
        cg = cigar.cigar_from_bam(h.cigar)
        #logger.debug(f"R: {rstart}-{rend}, len: {cg.get_len()}, Q: {qstart}-{qend}, len: {cg.get_len(ref=False)}")
        #print(h.mapq)
        if rstart > rend:
            logger.error(f"Inverted on Reference: {h}")
            continue
            # shouldn't ever occur, TODO maybe handle anyway?

        #logger.info(f"rstart: {rstart}, rend: {rend}, qstart: {qstart}, qend: {qend}")
        #logger.info(f"rtree: {str(reftree)}")
        #logger.info(f"qtree: {str(qrytree)}")

        rstartov = list(reftree[rstart])[0]
        qstartov = list(qrytree[qstart])[0]

        # simply append alignment if there is only one offset
        # as this happens quite often, this should save a lot of time
        # print(f'reftree: {reftree}, qrytree: {qrytree}, rend: {rend}, qend: {qend}')
        if rstartov == list(reftree[rend-1])[0] and qstartov == list(qrytree[qend-1])[0]:
            roff = rstartov.data
            qoff = qstartov.data
            aln = [rstart + roff, rend + roff, qstart + qoff, qend + qoff, rend - rstart +1, qend - qstart +1, cg.get_identity()*100, 1 if rstart < rend else -1, h.strand, rcid, qcid, cg.to_string()]

            # check to make sure alns match cigar length
            if aln[4] != cg.get_len(ref=True):
                    logger.error(f"Aln length not matching cigar length! {aln}")
            if aln[5] != cg.get_len(ref=False):
                    logger.error(f"Aln length not matching cigar length! {aln}")

            alns.append(aln)
            continue
        else:
            logger.warning(f"Multiple ({len(reftree[rstart:rend])}) offsets in one alignment!")

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
                aln = [rint.data + rstdel, rint.data + min(rend, rint.end) - rendel - max(rint.begin - rstart, 0),
                           qint.data + max(qstart, qint.begin), qint.data + min(qend, qint.end),
                           min(rend, rint.end) - rendel - rstdel - max(rint.begin - rstart, 0), min(qend, qint.end) - max(qstart, qint.begin),
                           qcg.get_identity()*100, 1 if rstart < rend else -1, 1 if qstart < qend else -1, rcid, qcid, qcg.to_string()]

                # check to make sure alns match cigar length
                if aln[4] != qcg.get_len(ref=True):
                    logger.error(f"Aln length not matching cigar length on ref! Occurred in {aln} (cg len {qcg.get_len(ref=True)})")
                if aln[5] != qcg.get_len(ref=False):
                    logger.error(f"Aln length not matching cigar length on qry! Occurred in {aln} (cg len {qcg.get_len(ref=False)})")
                alns.append(aln)

        # print('c')

    #logger.debug('Alignments traversed')
    # print('d')
    alns = pd.DataFrame(alns)
    #print(alns)
    if alns.empty:
        return None
    #print(alns[6])
    #alns[6] = alns[6].astype('float')

    alns = alns.loc[alns[6] > 90] # TODO: Alignment identity filter. This filter is not mandatory and the user might opt to remove this
    # count inverted alns as well
    #alns.loc[alns[8] == -1, 2] = alns.loc[alns[8] == -1, 2] + alns.loc[alns[8] == -1, 3]
    #alns.loc[alns[8] == -1, 3] = alns.loc[alns[8] == -1, 2] - alns.loc[alns[8] == -1, 3]
    #alns.loc[alns[8] == -1, 2] = alns.loc[alns[8] == -1, 2] - alns.loc[alns[8] == -1, 3]
    alns.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    alns.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    #print(alns[['aStart', 'aLen', 'bStart', 'bLen', 'iden']])
    # print('e')
    return None if alns.empty else alns

# </editor-fold>

cpdef generate_seqdict(fafin, mappingtrees, chrdict):
    return {org:('N'*_NULL_CNT).join([
        fafin[org].fetch(region = chrdict[org],
                         start = interval.data,
                         end = interval.data + interval.end - interval.begin).upper()
        for interval in sorted(mappingtrees[org])])
        for org in mappingtrees}


cpdef get_at_pos(alns, rchrom, rstart, rend, qchrom, qstart, qend):
    """
    Function that takes a Dataframe of alignments, and extracts the part of each alignment overlapping with the specified regions.
    :returns: Returns a new dataframe in the same format.
    """
    ret = deque()

    #logger.debug(f"printing ALNS matching this pos: {rchrom}, {rstart}-{rend}/{qstart}-{qend}")
    #print(alns.loc[(alns['achr'] == rchrom) & (alns['astart'] <= rstart) & (alns['aend'] >= rend) & (alns['adir'] == 1) &
    #                    (alns['bchr'] == qchrom) & (alns['bstart'] <= qstart) & (alns['bend'] >= qend) & (alns['bdir'] == 1)][['achr', 'bchr', 'astart', 'aend', 'bstart', 'bend', 'alen', 'blen']])


    # iterate over all alns overlapping on both qry and ref
    alns = get_overlapping(get_overlapping(alns, rstart, rend, chrom=rchrom), qstart, qend, chrom=qchrom, ref=False)
    #logger.debug(f"called with {rchrom}, {rstart}, {rend}, {qchrom}, {qstart}, {qend}, found overlapping {alns}")
    if alns is None:
        return None
    for _, aln in alns.iterrows():
        # cut off alignments to only in the gap we are realigning
        cg = cigar.cigar_from_string(aln.cg)

        # check that aln lengths are correct
        #if cg.get_len() != aln.aend - aln.astart + 1:
        #    logger.error(f"CIGAR len ({cg.get_len()}) not matching len on reference ({aln.aend - aln.astart + 1})!")
        #if cg.get_len(ref=False) != aln.bend - aln.bstart + 1:
        #    logger.error(f"CIGAR len ({cg.get_len()}) not matching len on reference ({aln.aend - aln.astart + 1})!")

        #logger.debug(f"Removing {rstart - aln.astart}, {aln.aend - rend} from aln with len {cg.get_len()}")
        #print(cg.to_string())
        srem, erem, cg = cg.trim(max(0, rstart - aln.astart), max(0, aln.aend - rend))
        
        # check that the positions after removing match
        #if srem != qstart - aln.bstart:
        #    logger.error(f"Mismatch during alignment trimming, start does not map on query! Should have removed {qstart - aln.bstart}, actually removed {srem}. CIGAR: {cg.to_string()}")
        #if erem != aln.bend - qend:
        #    logger.error(f"Mismatch during alignment trimming, end does not map on query! Should have removed {aln.bend - qend}, actually removed {erem}. CIGAR: {cg.to_string()}")

        ## check that lengths match
        #if rend - rstart + 1 != cg.get_len(ref=True):
        #    logger.error(f"Coordinate length ({rend - rstart + 1}) not matching cigar length ({cg.get_len(ref=True)}) on ref! Occurred in {aln}")
        #if qend - qstart + 1 != cg.get_len(ref=False):
        #    logger.error(f"Coordinate length ({qend - qstart + 1}) not matching cigar length ({cg.get_len(ref=False)}) on qry! Occurred in {aln}")

        # use cigar lens to force eager trimming of CIGARS
        # otherwise, I/D records at the end of the ALN could stick around, confusing later steps
        ret.append([rstart, rstart + cg.get_len(), qstart, qstart + cg.get_len(ref=False), cg.get_len(), cg.get_len(ref=False), cg.get_identity()*100,
                    aln.adir, aln.bdir, aln.achr, aln.bchr, cg.to_string()]) # we change neither orientation nor chromosome of the ALN


    if len(ret) == 0: # return no aln if none found
        return None

    return pd.DataFrame(ret, columns = ["astart", "aend", "bstart", "bend", "alen", "blen", "iden", "adir", "bdir", "achr", "bchr", 'cg'])

cdef syrify(alnsdf):
    if alnsdf is None:
        return None
    alnsdf.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    alnsdf.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    return alnsdf

cdef get_overlapping(alnsdf, start, end, chrom=None, ref=True, dir=1):
    """Helper Fn to filter an alignment DF for any region overlapping with [start:end] on a (default) or b (if `ref=False` is passed). If `dir` is passed and not 0 (default 1), also filter for the alignment direction (-1 = inverted, 1 non-inverted).
    Filters for a chromosome on a or b if specified, otherwise ignores chromosomes"""
    #return alnsdf.loc[!((alnsdf.bstart < start) ^ (alnsdf.bend > end))]
    if alnsdf is None or alnsdf.empty:
        return None
    startcol = alnsdf.astart if ref else alnsdf.bstart
    endcol = alnsdf.aend if ref else alnsdf.bend
    ret = alnsdf.loc[((alnsdf['achr' if ref else 'bchr'] == chrom) if chrom else True) &
                    ((alnsdf['adir' if ref else 'bdir'] == dir) if dir != 0 else True) & (
                    ((startcol >= start) & (endcol <= end)) | # get regions fully contained
                    ((startcol <= start) & (endcol >= end)) | # or starting before and ending beyond
                    ((startcol >= start) & (startcol < end)) | # or starting and ending beyond
                    ((endcol > start) & (endcol <= end) ) # or starting before and ending within
                    )]
    #print(f"Called with {start}, {end}, {chrom}, {ref}, {dir}")
    #print(f"{ret}")
    return ret

cpdef get_nonsyn_alns(alnsdf, reftree, qrytree):
    """
    Function that extracts alignments of sequence that has not been called as merasyntenic yet from a set of alignments, in preparation for the synteny identification part of realignment.
    This Fn assumes the input alignments are all on the same chromosome in the same direction and will report alignments corresponding to any position on the reference – these conditions are ensured by calling get_at_pos on alnsdf first. 
    :args:
    :alnsdf: Dataframe of alignments (eg produced by io.read_alnsfile).
    :reftree: An Intervaltree with a start coordinate for each region that has not been identified as merasyntenic yet in the chosen reference. Produced for all samples at once by construct_mts.
    :qrytree: An Intervaltree with a start coordinate for each region that has not been identified as merasyntenic yet in the query sequence. Produced for all samples at once by construct_mts.
    :returns: A Dataframe in the same format. If there are multiple non-adjacent non-merasyn segments in the tree, it may have more alignments than in the input, by splitting larger alns per region.
    """

    ret = []
    #logger.debug(f"get_nonsyn called with ref {reftree}, qry {qrytree}")
    for rint in reftree:
        # pre-fetch overlapping alns
        # do not drop now, another copy is probably slower anyway
        rintlen = rint.end - rint.begin
        rintalns = get_overlapping(alnsdf, rint.data, rint.data + rintlen)
        #rintalns = get_at_pos(alnsdf, None, rint.data, rint.data + rintlen, None, None)
        #logger.debug(f"{rint}: found {rintalns}")

        for qint in qrytree:
            qintlen = qint.end - qint.begin
            #logger.debug(str(get_at_pos(alnsdf, None, rint.data, rint.data + rintlen, None, qint.data, qint.data + qintlen)))
            ret.append(get_at_pos(rintalns, None, rint.data, rint.data + rintlen, None, qint.data, qint.data + qintlen))

    #logger.debug(f"Found: {ret}")
    if len(ret) == 0 or all([r is None for r in ret]):
        logger.warning(f"No alignments found in this region! This could be a repetitive region, or the alignments could be truncated!")
        return None
    return syrify(pd.concat(ret))


cpdef realign(df, qrynames, fastas, MIN_REALIGN_THRESH=None, MAX_REALIGN=None, NULL_CNT=None, mp_preset='asm10', ncores=1, pairwise=None, output_only_realign=False):
    """
    Function to find gaps between two coresyn regions and realign them to a new reference.
    Discovers all merasynteny.
    
    :arguments: A DataFrame with core and merasyn regions called by find_pansyn and the sample genomes and names. `mp_preset` designates which minimap2 alignment preset to use.
    :returns: A DataFrame with the added non-reference merasynteny
    """
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
    # init stuff
    ret = deque()#pd.DataFrame()
    n = len(qrynames) + 1 # to account for reference
    if not n == len(fastas) + 1:
        logger.error(f"More/less query names than fastas passed to process_gaps: {qrynames}, {fastas}")
        raise ValueError("Wrong number of fastas!")

    # load fasta files
    fafin = {qrynames[i]: pysam.FastaFile(fastas[i]) for i in range(len(qrynames))}

    # iterate through each gap between coresyn blocks
    # call the alignment/ functionality and merge   

    syniter = df.iterrows()
    try:
        syn = next(syniter)[1][0]
        # skip to first core
        while syn.get_degree() < n:
            if not output_only_realign:
                ret.append(syn)
            syn = next(syniter)[1][0]
            # print(vars(syn))
        old = syn
        # TODO: Misses the merasyn before the first coresyn region? This would become if there are no or very few coresyn regions
        # leon: yes, the merasyn before the first and after the last coresyn of a chromosome will be missed.
        # my thinking was that this would correspond to the highly polymorphic tails of the chromosomes
        # this seems to work out in Ampril (the first coresyn starts at <2 kb in my local test dataset),
        # but might cause problems in datasets with very little coresyn
        merasyns = dict()
        # CNT = 0
        while True:
            # CNT += 1
            # print(CNT)
            # if CNT == 100:
            #     break
            # find block between two coresyn regions
            merasyns = dict()

            # store merasyn-regions, to be added once we know if we need to realign
            refmerasyns = []
            while syn.get_degree() < n:
                refmerasyns.append(syn)
                syn = next(syniter)[1].iloc[0]

            merasyns[old.ref.org] = refmerasyns

            # syn must be core now

            # start and end of the non-ref region, on the reference
            end = syn.ref.start -1 # end/start are inclusive
            start = old.ref.end +1
            #logger.debug(f"Realigning between {start} and {end}. Borders on ref: {old.ref}, {syn.ref}")#\n Full borders {old}, {syn}")

            # if there is not enough novel sequence on any organism to realign, skip this realignment preemptively
            if end - start < _MIN_REALIGN_THRESH:
                if all([syn.ranges_dict[org].start - old.ranges_dict[org].end < _MIN_REALIGN_THRESH \
                        for org in syn.ranges_dict]):
                    if not output_only_realign:
                        ret.append(syn)
                    old = syn
                    syn = next(syniter)[1][0]
                    continue

            # between chromosomes, there isn't a single gap
            #print(f"Chrs: {syn.ref.chr}, {old.ref.chr}")
            if syn.ref.chr != old.ref.chr:
                logger.error("Chr case found: {syn.ref}, {old.ref}")
                old = syn
                if not output_only_realign:
                    ret.append(syn)
                syn = next(syniter)[1][0]
                continue

            #########
            ## Block has been extracted and is long enough;
            ## extract appropriate sequences, respecting already found merasyn
            #########

            #print(old, syn, syn.ref.start - old.ref.end)
            #print(merasyns)

            #TODO parallelise everything after this, enable nogil

            
            ## construct the mapping, and prepare sequences for realignment
            mappingtrees = construct_mts(refmerasyns, old, syn)
            seqdict = generate_seqdict(fafin, mappingtrees, {org: syn.ranges_dict[org].chr for org in mappingtrees})

            #if any([len(mappingtrees[org]) > 3 for org in mappingtrees]):
            #    logger.info(mappingtrees)
            #    logger.info(seqdict)

            if not seqdict: # if all sequences have been discarded, skip realignment
                #logger.info("Not aligning, not enough non-reference sequence found!")
                old = syn
                if not output_only_realign:
                    ret.append(syn)
                syn = next(syniter)[1][0]
                continue

            ## Realign iteratively until all synteny is found
            while len(seqdict) > 2: # realign until there is only one sequence left

                # TODO: Have some heuristic terminate realignment in highly repetitive regions
                # stop realignment if we have already found _MAX_REALIGN haplotypes
                if _MAX_REALIGN > 0 and len(merasyns) > _MAX_REALIGN:
                    break

                ## choose a reference
                # uses the sample containing the most non-merasyntenic sequence
                # if a dict of pairwise alns is passed, will always prefer samples in the dict
                if pairwise:
                    ref = max([(len(v) if k in pairwise else (-1)/len(v), k) for k,v in seqdict.items()])[1]
                else:
                    ref = max([(len(v), k) for k,v in seqdict.items()])[1]

                #print('ref:', ref)
                #print('On ref:', syn.ref.chr, start, end, end - start)
                #print({org:len(seq) for org, seq in seqdict.items()})
                #        print(org, ':', seqdict[org])

                refseq = seqdict[ref]
                del seqdict[ref]
                reftree = mappingtrees[ref]
                del mappingtrees[ref]

                refstart = old.ref.end if ref == 'ref' else old.ranges_dict[ref].end
                refend = syn.ref.start if ref == 'ref' else syn.ranges_dict[ref].start

                if True:#refend - refstart > 5000 and len(syns) < 4:
                    print(f"\n===\n>{ref} {list(reftree)}\n{refseq}")
                    print('\n'.join([f">{id} {list(mappingtrees[id])}\n{seq}" for id, seq in seqdict.items()]))



                logger.debug(f"Realigning {old.ref.chr}:{old.ref.end}-{syn.ref.start} (len {util.siprefix(syn.ref.start - old.ref.end)}) to {ref}. Seqdict lens: {[(k, len(v)) for k,v in seqdict.items()]}")

                if refstart > refend:
                    logger.error(f"{refstart} after {refend}! Seqdict {[(k, len(v)) for k,v in seqdict.items()]}")
                    #continue


                ## get alignments to reference construct alignment index from the reference
                # TODO: The alignment step is a major performance bottleneck, specially when aligning centromeric regions. If the expected memory load is not high, then we can easily parallelise align_concatseqs using multiprocessing.Pool. Here, I have implemented it hoping that it should not be a problem. If at some point, we observe that the memory footprint increases significantly, then we might need to revert it back.
                if pairwise and ref in pairwise:
                    # if we have pairwise alns, fetch them
                    logger.debug(f"Fetching from existing alignments. Left core: {old.ref} ({old.ranges_dict}). Right core: {syn.ref} ({syn.ranges_dict}). Ref {ref}")
                    
                    refalnsdict = pairwise[ref]
                    # get all the alns overlapping this region; syri should do the rest
                    # regions not in seqdict will be ignored
                    alns = {org: get_nonsyn_alns(
                                    get_at_pos(refalnsdict[org], old.ref.chr, refstart, refend, old.ranges_dict[org].chr, old.ranges_dict[org].end, syn.ranges_dict[org].start), # pre-process alignments to restrict to this realn region
                                    reftree, mappingtrees[org])
                            for org in seqdict}
                else:
                    # otherwise realign ourselves
                    logger.debug(f"Starting Alignment. Left core: {old.ref}. Right core: {syn.ref}. Ref {ref}")
                    if len(refseq) > 50000:
                        logger.debug(f"Starting parallel Alignment to {ref} between {refstart} and {refend} (len {util.siprefix(refend - refstart)})")
                        alignargs = [[seqdict[org], syn.ranges_dict[org].chr, mappingtrees[org]] for org in seqdict.keys()]
                        with Pool(processes=ncores) as pool:
                            # pool.starmap(partial(foo, d='x'), alignargs)
                            alns = pool.starmap(partial(align_concatseqs, refseq=refseq, preset=mp_preset, rcid=syn.ref.chr, reftree=reftree, aligner=None), alignargs)
                        alns = dict(zip(list(seqdict.keys()), alns))
                    else:
                        aligner = get_aligner(seq=refseq, preset=mp_preset)
                        alns = dict()
                        # print('seq')
                        for org, seq in seqdict.items():
                            if seq == '': # skip empty sequences
                                alns[org] = None
                                continue

                            # TODO: Currently (12.03.2024), this seems to be the most time-consuming step
                            logger.debug(f"Processing alignments for {org} to {ref}. Seq len {len(seq)}.")
                            alns[org] = align_concatseqs(seq, syn.ranges_dict[org].chr, mappingtrees[org], refseq, mp_preset, syn.ref.chr, reftree, aligner=aligner)


                # filter out alignments only containing inversions
                for org in alns:
                    if alns[org] is not None and all(alns[org].bDir == -1):
                        logger.warning(f"{org} in alns only contains inverted alignments: \n{alns[org]}")
                        alns[org] = None

                #logger.info(f"None/empty in Alignments: {[org for org in alns if alns[org] is None]}")
                #print(ref, refseq)
                #print(seqdict)

                ## run syri
                logger.debug("Running syri")
                syris = syri_get_syntenic(alns)
                # TODO MG: Replaced getsyriout with the synteny identification method from syri. Consider parallelizing this because for repetitive regions this loop would be expensive as well.

                ## Match ALNs, in preparation for merasyn identification
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
                    logger.info(f"No synteny to {ref} was found!")
                    continue

                # syns should be sorted
                # TODO: MG Question: what is the use of this function?
                # leon: I've found that sometimes, syri calls would overlap a bit,
                # which the pansyn identification doesn't like
                # this checks this isn't the case and corrects it in case it finds it
                # should probably not be strictly necessary here as we call syri ourselves
                # but doesn't hurt to check either

                
                # print(syns)
                ## Find merasyn in the syri calls
                pansyns = syns
                pansyns = pansyn.reduce_find_overlaps(syns, cores=1)

                # no need to recalculate the tree if no pansynteny was found
                if pansyns is None or pansyns.empty:
                    logger.info("No multisynteny was found in this round!")
                    continue

                # Add all merasyns with alphabetical sorting by reference name
                merasyns[ref] = [psyn[1][0] for psyn in pansyns.iterrows()]
                added = sum([len(x.ref) for x in merasyns[ref]])

                logger.info(f"Realigned {old.ref.chr}:{old.ref.end}-{syn.ref.start} (len {util.siprefix(syn.ref.start - old.ref.end)}) to {ref}. Found {util.siprefix(added)} (avg {util.siprefix(added/len(merasyns))}) of merasynteny.")

                ## recalculate mappingtrees from current merasyns to remove newly found merasynteny
                # TODO maybe directly remove, should be more efficient
                logger.debug(f"Old Mappingtrees: {mappingtrees}.\n Adding {merasyns[ref]}.")
                mappingtrees = subtract_mts(mappingtrees, merasyns[ref])
                logger.debug(f"New Mappingtrees: {mappingtrees}")
                # remove all orgs that have already been used as a reference
                for reforg in merasyns:
                    if reforg in mappingtrees:
                        del mappingtrees[reforg]

                ## extract the remaining sequences for future realignment

                # cache the current lens, for debugging
                lendict = {org: len(seqdict[org]) for org in seqdict}
                #logger.info("emit ending lens in seqdict")
                #logger.info(str(lendict))

                seqdict = generate_seqdict(fafin, mappingtrees, {org: syn.ranges_dict[org].chr for org in mappingtrees})
                if not seqdict: # if all sequences have been discarded, finish realignment
                    break

                # check that the sequence length has not been extended during the update
                # allow for up to one spacer to be inserted, though
                for org in seqdict:
                    if org in lendict: # eliminating sequences is always okay
                        #logger.info(f"Re-constructing {org} sequence. New len {util.siprefix(len(seqdict[org]))}, old {util.siprefix(lendict[org])}")
                        assert(len(seqdict[org]) <= lendict[org] + _NULL_CNT, "sequence length extended during update")


            # incorporate into output DF, sorted alphabetically by ref name
            # does nothing if no merasyn was found
            for org in sorted(merasyns.keys()):
                ret.extend(merasyns[org])

            # continue checking the next coresyn gap
            old = syn
            if not output_only_realign:
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

# Idea for additional fn
# finds private regions by scanning through the genome for regions not covered by any merasyn
# tracks current position along the genome
# challenge: non-coresyn regions
# => approach: sort merasyns by org, then subtract
# alternatively, use intervaltrees, subtract each merasyn
# maybe move to own file eventually?
# implement after refactoring of data structures
#cpdef find_private(syns, only_private=False):
#    # Finds 
#    pass

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

#
#cimport libc.stdio as cio
#cimport posix.unistd as unistd
#

## TODO: Make parameters adjustable. Also, now (12.03.24) this could be deprecated.
#cpdef getsyriout(coords, PR='', CWD='.', N=1, TD=500000, TDOLP=0.8, K=False, redir_stderr=False):
#    """DEPRECATED"""
#    BRT = 20
#    TUC = 1000
#    TUP = 0.5
#    T = 50
#    invgl = 1000000
#
#    #assert(len(list(np.unique(coords.aChr))) == 1)
#    try:
#        assert coords.aChr.nunique() == 1
#        assert coords.bChr.nunique() == 1
#    except AssertionError:
#        logger.error(f"Incorrect coords. More than one chromosome parsed. Ref chromosomes: {coords.aChr}. Qry chromosomes: {coords.bChr}")
#
#    cdef int oldstderr = -1
#    if redir_stderr:
#        #cio.fclose(cio.stderr)
#        #cio.stderr = cio.freopen(bytes(f"{CWD}/stderr", encoding='utf8'), "w", cio.stderr)
#        oldstderr = unistd.dup(unistd.STDERR_FILENO)
#        cio.freopen(bytes(f"{CWD}/stderr", encoding='utf8'), "w", cio.stderr)
#
#    # NOTE: syri requires that the coords table have same chromosome IDs for homologous chromosomes. When, the coords have different chromosome IDs, then manipulate the chroms IDs here
#    chromr = list(coords.aChr)[0] # there should only ever be one chr anyway
#    chromq = list(coords.bChr)[0]
#    samechrids = chromr == chromq
#    if not samechrids:
#        coords.bChr.replace(chromq, chromr, inplace=True)
#    chrom = chromr
#    # handle errors by return value; allows only showing output if there is a problem
#    # python errors coming after an error here will have normal stderr
#    # try:
#    #     # TODO: this function expects that the reference and query chromsome would have the same id. If that is not the case (pre-processing not done),then this function always crashes
#    #     syri(chrom, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD,tdolp=TDOLP)
#    # except ValueError:
#    #     print(coords[['aStart', 'aEnd', 'aLen', 'bStart', 'bEnd', 'bLen', 'iden', 'aDir', 'bDir']])
#    #     return None
#
#    if syri(chrom, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD, tdolp=TDOLP) == -1:
#        if redir_stderr:
#            logger.error("Redirecting stderr to console again")
#            #cio.fclose(cio.stderr)
#            #cio.stderr = oldstderr
#            unistd.close(unistd.STDERR_FILENO)
#            unistd.dup2(oldstderr, unistd.STDERR_FILENO)
#        logger.error("syri call failed on input:")
#        print(coords[['aStart', 'aEnd', 'aLen', 'bStart', 'bEnd', 'bLen', 'iden', 'aDir', 'bDir']])
#        if redir_stderr:
#            logger.error(f"syri stderr in '{CWD}/stderr'")
#        return None
#
#    #with multiprocessing.Pool(processes=N) as pool:
#    #    pool.map(partial(syri, threshold=T, coords=coords, cwdPath=CWD, bRT=BRT, prefix=PR, tUC=TUC, tUP=TUP, invgl=invgl, tdgl=TD,tdolp=TDOLP), chrs)
#
#    #TODO if runtime a problem: redo syri call to only call synteny => maybe configurable?
#    # Merge output of all chromosomes – still necessary for some reason
#    mergeOutputFiles([chrom], CWD, PR)
#
#    #TODO: Maybe not requires and can be removed?
#    # leon: outSyn fails if this isn't called, that's why it's left in there
#    # but yes, this step should be unnecessary
#    # In general, I think the syri calling should be done more elegantly --
#    # writing to and then reading from files is quite inefficient, especially
#    # for short realignments
#
#    #Identify meri-chromosomal events in all chromosomes simultaneously
#    getCTX(coords, CWD, [chrom], T, BRT, PR, TUC, TUP, N, TD, TDOLP)
#
#    # Recalculate syntenic blocks by considering the blocks introduced by CX events
#    outSyn(CWD, T, PR)
#    o = getsrtable(CWD, PR)
#    if not samechrids:
#        o.bchr.replace(chromr, chromq, inplace=True)
#
#    if redir_stderr:
#        #cio.fclose(cio.stderr)
#        #cio.stderr = oldstderr
#        unistd.close(unistd.STDERR_FILENO)
#        unistd.dup2(oldstderr, unistd.STDERR_FILENO)
#
#    if not K:
#        for fin in ["synOut.txt", "invOut.txt", "TLOut.txt", "invTLOut.txt", "dupOut.txt", "invDupOut.txt", "ctxOut.txt", "sv.txt", "notAligned.txt", "snps.txt"]:
#            try:
#                os.remove(CWD+PR+fin)
#            except OSError as e:
#                if e.errno != 2:    # 2 is the error number when no such file or directory is present https://docs.python.org/2/library/errno.html
#                    raise
#    return o
## END
#
