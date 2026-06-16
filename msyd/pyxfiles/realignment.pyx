#%%cython
#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import parasail
import logging
from collections import deque, defaultdict
from functools import partial
from multiprocessing import Pool

from intervaltree import IntervalTree, Interval

# I added these lines to hide all of the INFO logs from syri. If those are required then these lines can be removed
logging.getLogger('syri').setLevel(logging.WARNING)
logging.getLogger('getCTX').setLevel(logging.WARNING)

#from syri.synsearchFunctions import syri, mergeOutputFiles, outSyn, apply_TS, alignmentBlock, getSynPath
#from syri.tdfunc import getCTX
#from syri.writeout import getsrtable
from syri.synsearchFunctions import apply_TS, alignmentBlock, getSynPath

import msyd.util as util
import msyd.cigar as cigar
import msyd.intersection as intersection
import msyd.priv as priv
#import msyd.io
from msyd.multisyn import Multisyn, Private, MultisynContainer
from msyd.coords import Range


cdef:
    int _MIN_REALIGN_LEN = 100 # min length to realign regions
    int _MIN_SYN_ID = 80 # minimum % identity for a region to be considered syntenic
    int _MAX_REALIGN = 0 # max number of haplotypes to realign to; set to 0 to realign without limit
    int _MAX_HANGOVER = 2000 # max len of padding with neighbouring coresyn region to use during alignment
    int _SPACER_LEN = 100 # number of separators to use between blocks during alignment
    int _MIN_PRIV_THRESH = intersection.get_min_syn_thresh()
    int _GAP_OPEN = 6
    int _GAP_EXTEND = 2
    object _MATRIX = parasail.matrix_create("ACGTN", 2, -1)

# force Ns to never align
for i in range(5): # slice indexing does not work for matrices
    _MATRIX[4, i] = -100
    _MATRIX[i, 4] = -100

logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.DEBUG)

### Log total added bases
cdef volatile int ADDED_LEN = -1 # global counter, for logging; this isn't perfectly thread safe but that should be fine as we're just logging
# setting to -1 disables logging
cpdef int get_added_len():
    global ADDED_LEN
    return ADDED_LEN

cpdef start_log_added_len():
    global ADDED_LEN
    ADDED_LEN = 0

cpdef stop_log_added_len():
    global ADDED_LEN
    ADDED_LEN = -1

cpdef set_aln_params(int gap_open, int gap_extend, object matrix):
    global _GAP_OPEN, _GAP_EXTEND, _MATRIX
    _GAP_OPEN = gap_open
    _GAP_EXTEND = gap_extend
    _MATRIX = matrix

### Example
## AAAANNNNNBBBBBBB
## 4 bp seq, 5 bp spacer, 8 bp seq
## => [0, 3] -> A, [9, 15] -> B

cpdef nonsyns_to_mt(nonsyns):
    """
    Transforms a list of alignments for each organism into a dictionary of intervaltrees mapping an offset in the alignment reference to each position in the sample genome.
    """
    pos = 0
    tree = IntervalTree()
    for rng in nonsyns:
        if len(rng) > _MIN_REALIGN_LEN: # filter again for sufficient len
            tree[pos: pos + len(rng)] = rng.start # end is non-inclusive
            pos += len(rng) + _SPACER_LEN # add interval + spacer

    return tree

cpdef extract_nonsynsdict(merasyns, gap_intervals):#prevcore, nextcore):
    """
    Makes a dictionary containing an intervaltree with an offset mapping for each org containing enough non-aligned sequence to realign.
    Merasyns need to be sorted by position on reference.
    For each tree, the sequence in genome `org` at position `tree[pos].data - tree[pos].begin + pos` corresponds to the position `pos` in the synthetic query sequence.
    """
    listdict = defaultdict(list)
    offsetdict = {org:rng.start for org, rng in gap_intervals.items()} # stores the current offset in each org

    for merasyn in iter(merasyns):
        # mark private regions as covered, if they are 
        if merasyn.is_private():
            if merasyn.ref.org != 'ref':
                offsetdict[merasyn.ref.org] = merasyn.ref.end +1
            continue

        # iterate through all multisyns found so far
        for org, rng in merasyn.ranges_dict.items():
            if not org in offsetdict: # no need to process
                continue
            #print(f"{offsetdict[org]}, {posdict[org]}, {rng}, {mtrees[org]}")
            l = rng.start - offsetdict[org] # len of the region to be added.
            if l < 0:
                logger.error(f"improper sorting: {rng.start} < {offsetdict[org]}") # improper sorting – skip
                continue

            # check if this interval would be redundant
            if org in listdict:
                prev = listdict[org][-1]
                #print(prev[0] + prev[1], offsetdict[org])
                if prev.end == offsetdict[org]: 
                    # check if offset + len matches the current offset; extend prev interval instead
                    # no +1 because the end is not inclusive
                    prev.end = prev.end + l
                elif -10 < prev.end - offsetdict[org] < 10: # log to detect off by ones
                    logger.debug(f"Close miss in {prev}, {offsetdict[org]}")

            if l > _MIN_REALIGN_LEN: # otherwise add to the tree if it's large enough
                listdict[org].append( Range(org, gap_intervals[org].chrom, offsetdict[org], offsetdict[org] + l - 1) )

            # all up to the end of this region has been added
            offsetdict[org] = rng.end + 1

    # see if there's any sequence left to realign after processing the merasyn regions
    for org, offset in offsetdict.items():
        l = gap_intervals[org].end - offset
        if l >= _MIN_REALIGN_LEN:
            listdict[org].append( Range(org, gap_intervals[org].chrom, offsetdict[org], offsetdict[org] + l ) )

    return listdict #listdict_to_mts(listdict)
# END

cpdef subtract_nonsynsdict(nonsyns_dict, merasyns):
    # merasyns has to be sorted
    cdef:
        cur_index = defaultdict(int)
        cur_pos = defaultdict(int)
        ret = defaultdict(list)

    # subtract merasyns
    for msyn in iter(merasyns):
        for org, rng in msyn.iter_orgs_ranges():
            # skip until msyn is covered
            org_ind = cur_index[org]
            org_nonsyns = nonsyns_dict[org]
            while org_ind < len(org_nonsyns):# and org_nonsyns[org_ind].end < rng.end:
                org_rng = org_nonsyns[org_ind]
                if org_rng.end < rng.start: # fully before, nonoverlapping
                    if org_rng.end > cur_pos[org]: # skip if fully within curpos
                        ret_rng = Range(org_rng.org, org_rng.chrom, max(org_rng.start, cur_pos[org]), org_rng.end)
                        if len(ret_rng) > _MIN_REALIGN_LEN:
                            ret[org].append(ret_rng)
                # left overlap
                elif org_rng.start < rng.start: # left overlap
                    leftov = Range(org_rng.org, org_rng.chrom, max(cur_pos[org], org_rng.start), rng.start - 1)
                    if len(leftov) > _MIN_REALIGN_LEN:
                        ret[org].append(leftov)
                    cur_pos[org] = rng.end + 1 # 
                #elif org_rng.start > rng.start and org_rng.end < rng.end: # fully contained, do not add
                #    pass
                # break if the merasyn ends, otherwise increment index
                if org_rng.end > rng.end:
                    break
                else:
                    orgind += 1
            # if ended due to lack of index, no further region to add
            # update counter
            cur_index[org] = org_ind
            cur_pos[org] = rng.end + 1
    return ret


cdef compute_gaps(chrom: str, prevcore: Multisyn, nextcore: Multisyn, lendict: dict):
    """
    Extracts the length of the gap on each organism.
    Can handle the start and end case where prev/nextcore are None.
    Takes the maximal chromosome length info and orgs from lendict, which can be obtained by callingget_len_dict on the sequence handler.
    """
    cdef ret = dict()

    for org in lendict: # prevcore and nextcore may both be None if at start/end
        start = prevcore.get_range(org).end +1 if prevcore else 0
        end = nextcore.get_range(org).start -1 if nextcore else lendict[org]
        assert end - start >= -1 # -1 is a gap of 0

        #if end - start > _MIN_REALIGN_LEN:
        ret[org] = Range(org=org, chrom=chrom, start=start, end=end)

    return ret


cpdef realign(chrcont, qrynames, seqh, MIN_REALIGN_LEN=None, MIN_SYN_ID=None, MAX_REALIGN=None, SPACER_LEN=None, mp_preset='asm20', ncores=1, annotate_private=True, pairwise=None, output_only_realign=False, debug_export=False):
    """
    High-level interface to the realignment functionality.
    Takes a dict of DataFrames containing per-chromosome Multisyn annotations.
    Gaps between core syntenic regions are iteratively realigned to identify non-reference merasynteny.

    :params:
    :df: DataFrame containing existing Multisyn annotations
    :qrynames: Names of the non-ref. sequences
    :fastas: Dictionary mapping the name of a sample to its genome as a pysam FastaFile object
    :mp_preset: minimap2 alignment preset to use. Default `asm20`.
    :annotate_private: After a sequence has been used as a reference in realignment and no synteny has been found, we can say that it is structurally private. If this flag is set, it will be annotated as such. Default `True`.
    :ncores: Number of cores to use during the realignment. Default `1`.
    :pairwise: Optionally, a dictionary of dictionaries containing pysam AlignmentFile objects for all pairwise alignments. Useful if your sequences have already been all vs all aligned, or you want to use a different aligner than minimap2.
    :output_only_realign: Flag on whether to also output annotations not obtained through realignment. Default `False`. Mostly intended for debugging.

    :returns: A DataFrame of Multisyn objects corresponding to the annotations after realignment.
    """
    # export globals if changes to configs are specified
    if MIN_REALIGN_LEN is not None and MIN_REALIGN_LEN >= 0:
        global _MIN_REALIGN_LEN
        _MIN_REALIGN_LEN = int(MIN_REALIGN_LEN)
    if MIN_SYN_ID is not None and MIN_SYN_ID >= 0:
        global _MIN_SYN_ID
        _MIN_SYN_ID = int(MIN_SYN_ID)
    if MAX_REALIGN is not None and MAX_REALIGN >= 0:
        global _MAX_REALIGN
        _MAX_REALIGN = int(MAX_REALIGN)
    if SPACER_LEN is not None and SPACER_LEN >= 0:
        global _SPACER_LEN
        _SPACER_LEN = int(SPACER_LEN)

    start_log_added_len()

    _process_gaps = partial(process_gaps, seqh=seqh, annotate_private=annotate_private, pairwise=pairwise, output_only_realign=output_only_realign, debug_export=debug_export)
    ret = chrcont.apply_chroms_par(_process_gaps, ncores=ncores)

    logger.info(f"Realigment done. Added {util.siprefix(get_added_len())} in total.")
    stop_log_added_len()

    return ret

#NOTE cpdef'd to enable using functools.partial.
#NOTE Consider wrapping with cython to enable re-cdefing this?
cpdef process_gaps(chrom:str, msyncont:MultisynContainer, seqh:seq.SeqHandler, alnparams=None, ncores=1, annotate_private=True, pairwise=None, output_only_realign=False, debug_export=False):
    """
    Workhorse function of the realignment functionality.
    Takes a DF of multisyns, finds gaps of sufficient size between coresyn regions in the DF to process.
    In these gaps, existing merasynteny is masked, and a new reference is chosen.
    All sequences that have not previously been used as reference are realigned to the new reference, and pairwise synteny is called.
    On the pairwise synteny to the new reference, the synteny intersection algorithm is run, identifying new merasynteny.
    This is iterated until either there are no regions without annotated merasynteny annotation left, but at most `MAX_REALIGN` times.

    :params:
    :df: DataFrame containing existing Multisyn annotations
    :qrynames: Names of the non-ref. sequences
    :fastas: Dictionary mapping the name of a sample to its genome as a pysam FastaFile object
    :mp_preset: minimap2 alignment preset to use. Default `asm20`.
    :annotate_private: After a sequence has been used as a reference in realignment and no synteny has been found, we can say that it is structurally private. If this flag is set, it will be annotated as such. Default `True`.
    :ncores: Number of cores to use during the realignment. Default `1`.
    :pairwise: Optionally, a dictionary of dictionaries containing pysam AlignmentFile objects for all pairwise alignments. Useful if your sequences have already been all vs all aligned, or you want to use a different aligner than minimap2.
    :output_only_realign: Flag on whether to also output annotations not obtained through realignment. Default `False`. Mostly intended for debugging.

    :returns: A DataFrame of Multisyn objects corresponding to the annotations after realignment.
    """
    #print(chrom, msyncont, seqh)
    # init stuff
    cdef:
        list ret = list()#deque()#pd.DataFrame()
        orgs = msyncont.orgs
        dict lendict = seqh.get_len_dict(chrom)

    if alnparams:
        #alnparams* apparently not working
        set_aln_params(alnparams[0], alnparams[1], alnparams[2])

    # iterate through each gap between coresyn blocks
    # call the alignment/ functionality and merge   
    prevcore = None # store a backlink to the previous coresyn
    for nextcore, merasyns in msyncont.iter_cores_acc():
        gap_intervals = compute_gaps(chrom, prevcore, nextcore, lendict)
        logger.info(f"Found gap {gap_intervals}")
        # debugging, skip large segments
        gap_intervals = {org:gap for org, gap in gap_intervals.items() if _MIN_REALIGN_LEN <= len(gap) <= 20000}
        #logger.warning(f"Skipping {gap} on {org} (too long)!")

        # Realign the gap, if it has a region larger than _MIN_REALIGN_LENGTH
        if gap_intervals:
            logger.info(f"Realigning after filtering: {gap_intervals}")
            ##NOTE implement deduplication here
            # make a dict of representative seqs + offset from gap_intervals?

            ## Find nonsyns, iterate realignment
            nonsyns_dict = extract_nonsynsdict(merasyns, gap_intervals)
            realsyns = iterate_reprocessing(nonsyns_dict, seqh, ncores=ncores, pairwise=pairwise, annotate_private=annotate_private)

            ##NOTE reduplicate afterwards

            # write directly if only storing realigned;
            # otherwise insert into DF respecting sorting
            if output_only_realign:
                ret.extend(sorted(realsyns))
            else:
                merasyns.extend(realsyns)

        # DONE with realignment

        logger.debug(f"{merasyns}")
        # add multisyns if 
        if not output_only_realign:
            if prevcore:
                ret.append(prevcore)
            ret.extend(sorted(merasyns))
        prevcore = nextcore

    # Done with all gaps, return as new MultisynContainer
    return MultisynContainer.from_iterable(ret)
# END

cdef iterate_reprocessing(nonsyns_dict, seqh, ncores=1, pairwise=None, annotate_private=False):
    global ADDED_LEN
    ## construct the mapping, and prepare sequences for realignment
    cdef:
        list ret = list()
        list added_lens = list()
        list added_privs = list()
        list used_refs = list()

    ## Realign iteratively until all synteny is found
    while True:
        # fetch sequences
        if not nonsyns_dict: # if all remaining are too small
            break
        logger.debug(f"Nonsyns: {list(nonsyns_dict.items())}")

        ## choose a reference
        # uses the sample containing the most non-synteny
        # if pairwise is passed, prefer samples in the dict
        if pairwise:
            ref = max([(sum([len(nonsyn) for nonsyn in nonsyns]) if org in pairwise else (-1)/sum([len(nonsyn) for nonsyn in nonsyns]), org) for org, nonsyns in nonsyns_dict.items()])[1]
        else:
            ref = max([(sum([len(nonsyn) for nonsyn in nonsyns]), org) for org, nonsyns in nonsyns_dict.items()])[1]

        ## assemble reference concatseq & mappingtree
        ref_concatseq = ('N'*_SPACER_LEN).join([seqh.get_range(rng) for rng in nonsyns_dict[ref]])
        ref_mt = nonsyns_to_mt(nonsyns_dict[ref])#TODO

        ## align orgs
        # align nonsyns to ref concatseq
        if pairwise:
            raise NotImplemented("pairwise aln support is not available currently.")
            #TODO

        syns_dict = dict()
        for org, nonsyns in nonsyns_dict.items(): #NOTE parallelize?
            alns = aln_nonsyns(ref_concatseq, ref_mt, nonsyns, seqh, ref, nonsyns_dict[ref][0].chrom)
            #alns = ret[0]
            #unalns = ret[1]

            logger.debug(f"{org}, aln: {alns}")#, unaln: {unalns}")
            if not alns:
                continue
            syns = syri_get_syntenic(ref, syrify(alns))
            if syns: # do not add empty ones
                syns_dict[org] = syns
        logger.debug(f"Synsdict: {list(syns_dict.items())}")

        # Find merasyn in the realignment syri calls
        msyns = None
        if syns_dict:
            msyns = intersection.reduce_find_overlaps(syns_dict.values(), cores=1)#ncores)
        logger.debug(f"Msyns: {msyns}")

        # recompute nonsyn regions, account for new multisynteny
        #TODO this can be done efficiently by subtracting from the old nonsyndict
        if msyns:
            nonsyndict = subtract_nonsynsdict(nonsyns_dict, msyns)

        """
        # hangovers were added, remove the left and right coresyn
        if _MAX_HANGOVER > 0: 
            if msyns and len(msyns) >= 2:
                msyns = MultisynContainer(init=msyns._backing[1:-1], orgs=msyns.orgs)
            else:
                logger.warning("Error in realignment (hangovers added, but unaligned)!")
                logger.info("This may indicate a highly repetitive region, or incorrect alignment parameters.")

            #assert len(msyns) >= 2
            #assert msyns[0].get_degree() == len(seqdict)
            #assert msyns[-1].get_degree() == len(seqdict)
            #msyns = MultisynContainer(init=msyns._backing[1:-1], orgs=msyns.orgs)
        """

        if annotate_private:
            # after aligning all against ref, we can call the remainder as private to ref
            added_privs.append(sum([len(x) for x in nonsyns_dict[ref]]))
            ret.extend([Private(rng) for rng in nonsyns_dict[ref]])
        # no more to discover on ref
        used_refs.append(ref)
        del nonsyns_dict[ref]

        ## log length of sequences, append to ret
        if msyns:
            added_lens.append(sum([len(x.ref) for x in iter(msyns)]))
            ret.extend(msyns)

        # counts all sequences that are still above _MIN_REALIGN_LENGTH
        if len(nonsyns_dict) <= 1 or len(used_refs) >= _MAX_REALIGN:
            break

    logger.info(f"Realigned gap_intervals. Found {[util.siprefix(a) for a in added_lens]} aligning to {used_refs}")
    # log globally how much sequence was found during realignment
    if ADDED_LEN >= 0:
        ADDED_LEN += sum(added_lens)

    if annotate_private:
        logger.info(f"Found {[util.siprefix(a) for a in added_privs]} of private sequence.")

    return ret

#    ## get alignments to reference construct alignment index from the reference
#    # if we have pairwise alns, fetch & prepare them
#    if pairwise and ref in pairwise:
#        logger.debug(f"Fetching {refrng} from existing alignments.")
#        
#        refalnsdict = pairwise[ref]
#        # get all the alns overlapping this region; syri should do the rest
#        # regions not in seqdict will be ignored
#        alns = {org: get_nonsyn_alns(
#                        get_at_pos(refalnsdict[org], refrng, gap_intervals[org]), # pre-process alignments to restrict to this realn region
#                        refmtree, mtrees[org])
#                for org in seqdict if not org == ref}
#       return alns

cpdef aln_nonsyns(str refconcat, object refmt, list nonsyns, object seqh, str refchrom, str reforg):
    """
    Aligns a list of nonsyntenic ranges to a concatenated reference sequence.
    Re-maps the positions to reference/organism space.
    :returns: a DF containing the alns
    """
    #NOTE allow multiple alns per nonsyn region?
    # maybe initially no, add later?
    #NOTE could switch over to avoid refconcat as well
    # => lower RAM requirements, probably faster
    # => report all alns, let syri choose
    cdef:
        list alns = list()
        list unaligned = list()

    #NOTE could parallelise on this level as well?
    for nonsyn in nonsyns:
        ## get sequence
        seq = seqh.get_range(nonsyn)
        #logger.debug(f"{seq[:100]}, {refconcat[:100]}")

        ## aln to ref concatseq
        # do a semiglobal alignment to allow matching the right ID on ref
        aln = None
        #NOTE refactor out into separate fn, impl dispatch to minimap2 for larg regions?
        # fn should return two ranges + cigar directly
        try:
            aln = parasail.sg_dx_trace_striped_32(seq, refconcat, _GAP_OPEN, _GAP_EXTEND, _MATRIX)
            #logger.debug(f"{aln}")
            if not aln or aln.score <= 0: # no alignment found
                unaligned.append(nonsyn)
                continue

            cg = cigar.cigar_from_string(str(aln.cigar.decode)) #NOTE if slow, use bytes directly
            iden = cg.get_identity() # floating point no
            if iden*100 < _MIN_SYN_ID: # no w/ high identity found
                #NOTE do full local alignment? split region?
                unaligned.append(nonsyn)
                continue
        except ValueError as ve:
            logger.warning(f"Error during parasail call: {ve}")
            #TODO this isn't skipping?
            # try catch entire block?
            unaligned.append(nonsyn)
            continue
        
        ## decode aln, map to appropriate pos

        # shouldn't be necessary for a semiglobal aln,
        # still cleaner though in case of switching to full local
        aln_qury = Range(nonsyn.org, nonsyn.chrom,
                         nonsyn.start + aln.cigar.beg_query,
                         nonsyn.start + aln.end_query)

        # remap reference pos
        startint = list(refmt[aln.cigar.beg_ref])[0]
        endint = list(refmt[aln.end_ref])[0]
        aln_ref = Range(reforg, refchrom,
                        startint.data + aln.cigar.beg_ref - startint.begin,
                        endint.data + aln.end_ref - endint.begin)

        # add aln
        logger.debug(f"{aln_ref}, {aln_qury}, {cg}")
        alns.append((aln_ref, aln_qury, cg))


        #TODO additional local aln step?
        # probably best to test if necessary first
    return alns#, unaligned

cdef syri_get_syntenic(reforg, alns):
    # Synteny call parameters
    # all except T are unused
    #BRT = 20
    #TUC = 1000
    #TUP = 0.5
    #invgl = 1000000
    cdef:
        int T = 50
        dict syns = {}

    # check Chrs
    if not alns.aChr.nunique() == 1 and alns.bChr.nunique() == 1:
        logger.error(
            f"Incorrect coords. More than one chromosome parsed. Ref chromosomes: {alns.aChr}. Qry chromosomes: {alns.bChr}")
        return None

    # NOTE: syri requires that the coords table have same chromosome IDs for homologous chromosomes. When, the coords have different chromosome IDs, then manipulate the chroms IDs here
    chromr = list(alns.aChr)[0]  # there should only ever be one chrom anyway
    chromq = list(alns.bChr)[0]
    samechrids = chromr == chromq
    if not samechrids:
        alns.bChr.replace(chromq, chromr, inplace=True)
    chromo = chromr
    logger.debug(f"alns: {alns}")

    coordsData = alns[(alns.aChr == chromo) & (alns.bChr == chromo) & (alns.bDir == 1)]
    syndf = apply_TS(coordsData.aStart.values, coordsData.aEnd.values, coordsData.bStart.values,
                  coordsData.bEnd.values, T)

    # clean up graph
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
            
    if (not syndf) or (not blocks):
        logger.info(f"All alignments filtered out!")
        # all alns filtered out
        return None

    # get path through the graph
    synPath = getSynPath(blocks)
    synData = coordsData.iloc[synPath].copy()

    if not samechrids:
        synData.bChr.replace(chromr, chromq, inplace=True)

    synData.columns = list(map(str.lower, synData.columns))

    # return early if there is no large-scale synteny
    MIN_SYN_THRESH = intersection.get_min_syn_thresh()
    if synData.empty or\
            (synData['aend'] - synData['astart']).sum() < MIN_SYN_THRESH or\
            (synData['bend'] - synData['bstart']).sum() < MIN_SYN_THRESH:
        logger.warning(f"No synteny found in realignment syri call!")
        return None

    # subset to only relevant columns for the realignment
    synData = synData[['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend', 'cigar']]

    # make into multisyn objects, store in dataframe
    #NOTE necessary to convert back to cigar?
    return MultisynContainer.from_iterable(
            [Multisyn(ref=Range(reforg, syn['achr'], syn['astart'], syn['aend']),
                     ranges_dict={org:Range(org, syn['bchr'], syn['bstart'], syn['bend'])},
                     cigars_dict={org:cigar.cigar_from_string(syn['cigar'])})
           for _, syn in synData.iterrows()])


#############################
#### Use all vs all alns ####
#############################
cpdef get_at_pos(alns, rrng, qrng):
    """
    Function that takes a Dataframe of alignments, and extracts the part of each alignment overlapping with the specified regions.
    :returns: Returns a new dataframe in the same format.
    """
    cdef:
        list ret = list()

    # iterate over all alns overlapping on both qry and ref
    alns = get_overlapping(get_overlapping(alns, rrng.start, rrng.end, chrom=rrng.chrom), qrng.start, qrng.end, chrom=qrng.chrom, ref=False)

    if alns is None:
        return None
    for _, aln in alns.iterrows():
        cg = cigar.cigar_from_string(aln.cg)

        # check that aln lengths are correct
        if cg.get_len() != aln.aend - aln.astart + 1:
            logger.error(f"CIGAR len ({cg.get_len()}) not matching len on reference ({aln.aend - aln.astart + 1})!")
        if cg.get_len(ref=False) != aln.bend - aln.bstart + 1:
            logger.error(f"CIGAR len ({cg.get_len()}) not matching len on reference ({aln.aend - aln.astart + 1})!")

        # trim alns to only the gap we are realigning, to make subsequent drops more efficient
        srem, erem, cg = cg.trim(max(0, rrng.start - aln.astart), max(0, aln.aend - rrng.end))
        #_, _, cg = cg.trim(max(0, rstart - aln.astart), max(0, aln.aend - rend))
        
        # check that the positions after removing match
        if srem != qrng.start - aln.bstart:
            logger.error(f"Mismatch during alignment trimming, start does not map on query! Should have removed {qrng.start - aln.bstart}, actually removed {srem}. CIGAR: {cg.to_string()}")
        if erem != aln.bend - qrng.end:
            logger.error(f"Mismatch during alignment trimming, end does not map on query! Should have removed {aln.bend - qrng.end}, actually removed {erem}. CIGAR: {cg.to_string()}")

        ## check that lengths match
        if len(rrng) != cg.get_len(ref=True):
            logger.error(f"Coordinate length ({len(rrng)}) not matching cigar length ({cg.get_len(ref=True)}) on ref! Occurred in {aln}")
        if len(qrng) != cg.get_len(ref=False):
            logger.error(f"Coordinate length ({len(qrng)}) not matching cigar length ({cg.get_len(ref=False)}) on qry! Occurred in {aln}")

        # use cigar lens to force eager trimming of CIGARS
        # otherwise, I/D records at the end of the ALN could stick around, confusing later steps
        ret.append([rrng.start, rrng.start + cg.get_len(), qrng.start, qrng.start + cg.get_len(ref=False), cg.get_len(), cg.get_len(ref=False), cg.get_identity()*100,
                    aln.adir, aln.bdir, aln.achr, aln.bchr, cg.to_string()]) # we change neither orientation nor chromosome of the ALN


    if len(ret) == 0: # return no aln if none found
        return None

    return pd.DataFrame(ret, columns = ["astart", "aend", "bstart", "bend", "alen", "blen", "iden", "adir", "bdir", "achr", "bchr", 'cg'])

cdef syrify(alns):
    """
    Helper fn to format alignment dfs into the format SyRI uses.
    """
    if not alns: # keep empty stuff empty
        return None
    alnsdf = pd.DataFrame([[refrng.start, refrng.end, qryrng.start, qryrng.end, len(refrng), len(qryrng), cg.get_identity()*100, 1, 1, refrng.chrom, qryrng.chrom, cg] for refrng, qryrng, cg in alns])
    #alnsdf = pd.concat([[aln[0].start, aln[0].end, aln[1].start, aln[1].end, len(aln[0]), len(aln[1]), aln[2].get_identity(), 1, 1, aln[0].chrom, aln[1].chrom, aln[2]] for aln in alns])
    alnsdf.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    alnsdf.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    logger.debug(f"Alnsdf: {alnsdf}")
    print(alnsdf.loc[0])
    return alnsdf

cdef syrify_df(alnsdf):
    """
    Helper fn to format alignment dfs into the format SyRI uses.
    """
    if alnsdf is None:
        return None
    alnsdf.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
    alnsdf.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    return alnsdf

cdef get_overlapping(alnsdf, start, end, chrom=None, ref=True, dir=1):
    """
    Helper Fn to filter an alignment DF for any region overlapping with [start:end] on a (default) or b (if `ref=False` is passed). If `dir` is passed and not 0 (default 1), also filter for the alignment direction (-1 = inverted, 1 non-inverted).
    Filters for a chromosome on a or b if specified, otherwise ignores chromosomes.
    """
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
    return ret

cpdef get_nonsyn_alns(alnsdf, reftree, qrytree):
    """
    Function that extracts alignments of sequence that has not been called as merasyntenic yet from a set of alignments, in preparation for the synteny identification part of realignment.
    This Fn assumes the input alignments are all on the same chromosome in the same direction and will report alignments corresponding to any position on the reference – these conditions are ensured by calling get_at_pos on alnsdf first. 
    :args:
    :alnsdf: Dataframe of alignments (eg produced by msyd.io.read_alnsfile).
    :reftree: An Intervaltree with a start coordinate for each region that has not been identified as merasyntenic yet in the chosen reference. Produced for all samples at once by construct_mts.
    :qrytree: An Intervaltree with a start coordinate for each region that has not been identified as merasyntenic yet in the query sequence. Produced for all samples at once by construct_mts.
    :returns: A Dataframe in the same format. If there are multiple non-adjacent non-merasyn segments in the tree, it may have more alignments than in the input, by splitting larger alns per region.
    """

    ret = []
    for rint in reftree:
        # pre-fetch overlapping alns
        # do not drop now, another copy is probably slower anyway
        rintlen = rint.end - rint.begin
        rintalns = get_overlapping(alnsdf, rint.data, rint.data + rintlen)

        for qint in qrytree:
            qintlen = qint.end - qint.begin
            ret.append(get_at_pos(rintalns, Range(start=rint.data, end=rint.data + rintlen), Range(start=qint.data, end=qint.data + qintlen)))

    #logger.debug(f"Found: {ret}")
    if len(ret) == 0 or all([r is None for r in ret]):
        logger.warning("No alignments found in this region! This could be a repetitive region, or the alignments could be truncated!")
        return None
    return syrify_df(pd.concat(ret))

## DEPRECATED
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
#    chromr = list(coords.aChr)[0] # there should only ever be one chrom anyway
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
#cdef get_aligner(seq, preset, ns=True):
#    #aligner = mp.Aligner(seq=refseq, preset=preset, scoring=[1, 19, 39, 81, 3, 1, 10])
#    # using values from https://github.com/lh3/minimap2/blob/9b0ff2418c298f5b5d0df12b137896e5c3fb0ef4/options.c#L134
#    # https://github.com/lh3/minimap2/issues/155
#    # https://github.com/lh3/minimap2/blob/0cc3cdca27f050fb80a19c90d25ecc6ab0b0907b/python/README.rst?plain=1#L93
#    # values from the manpage, under presets -> asm5
#    #-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50
#    
#    #aligner = mp.Aligner(seq=seq, preset=preset, scoring=[1, 19, 39, 81, 39, 81, 100]) if ns else mp.Aligner(seq=seq, preset=preset)
#    # set --score-N parameter to 10
#    aligner = mp.Aligner(seq=seq, preset=preset, sc_ambi=10, max_chain_skip=255)
#    return aligner
#
#cpdef aln_concatseqs(seq, refseq, preset, aligner=None):
#    """
#    Function to align the concatenated sequences as they are and then remap the positions to the positions in the actual genome.
#    Both sequences should be on the same chromosomes.
#    Splits alignments that span multiple offsets into one alignment per offset
#    """
#    # Parse aligner from parent function when not using multiprocessing.Pool. When using Pool, define aligner here
#    if aligner is None:
#        aligner = get_aligner(refseq, preset)
#
#    return aligner.map(seq, extra_flags=0x4000000) # this is the --eqx flag, causing X/= to be added instead of M tags to the CIGAR string
#
#cpdef translate_concatalns(matches, qcid, qrytree, rcid, reftree):
#    # traverse alignments
#    alns = deque()
#    #logger.debug(f"{list(m)}")
#    list(matches)
#    logger.debug("Raw alns:")
#    logger.debug(f"{list(matches)}")
#    #NOTE simplify if/when removing spacers
#    for h in matches:
#        rstart: int = h.r_st
#        rend: int = h.r_en -1 # use inclusive indices
#        qstart: int = h.q_st
#        qend: int = h.q_en -1 # use inclusive indices
#        cg = cigar.cigar_from_bam(h.cigar)
#
#        if rstart > rend:
#            # shouldn't ever occur, TODO maybe handle anyway?
#            logger.error(f"Inverted on Reference: {h}")
#            continue
#
#        rstartov = list(reftree[rstart])[0]
#        qstartov = list(qrytree[qstart])[0]
#
#        # shortcut to simply append alignment if there is only one offset
#        # as this happens quite often, this should save a lot of time
#        #print(reftree, rend, qrytree, qend)
#        if rstartov == list(reftree[rend-1])[0] and qstartov == list(qrytree[qend-1])[0]:
#            roff = rstartov.data
#            qoff = qstartov.data
#            aln = [rstart + roff, rend + roff, qstart + qoff, qend + qoff, rend - rstart +1, qend - qstart +1, cg.get_identity()*100, 1 if rstart < rend else -1, h.strand, rcid, qcid, cg.to_string()]
#
#            # check to make sure alns match cigar length
#            if aln[4] != cg.get_len(ref=True):
#                    logger.error(f"[simple case] Aln length {aln[4]} not matching cigar length on ref {cg.get_len(ref=True)}! Occurred in {aln}")
#            if aln[5] != cg.get_len(ref=False):
#                    logger.error(f"[simple case] Aln length {aln[5]} not matching cigar length on qry {cg.get_len(ref=False)}! Occurred in {aln}")
#
#            alns.append(aln)
#            continue
#
#        #NOTE simplify if/when removing spacers
#        # multiple offsets in alignment; split for each offset
#        for rint in sorted(reftree[rstart:rend]):
#            # subset alignment to this reference offset interval
#            qstdel, rcg = cg.get_removed(max(rint.begin - rstart, 0))
#            qendel, rcg = rcg.get_removed(max(rend - rint.end, 0), start=False)
#            for qint in sorted(qrytree[qstart + qstdel:qend - qendel]):
#                # subset to the query offset, respecting the subsetting done so far
#                rstdel, qcg = rcg.get_removed(max(qint.begin - qstdel - qstart, 0), ref=False)
#                rendel, qcg = qcg.get_removed(max(qend - qint.end - qendel, 0), ref=False, start=False)
#
#                aln = [rint.data + rstdel, rint.data + min(rend, rint.end) - rendel - max(rint.begin - rstart, 0),
#                           qint.data + max(qstart, qint.begin), qint.data + min(qend, qint.end),
#                           min(rend, rint.end) - rendel - rstdel - max(rint.begin - rstart, 0), min(qend, qint.end) - max(qstart, qint.begin),
#                           qcg.get_identity()*100, 1 if rstart < rend else -1, 1 if qstart < qend else -1, rcid, qcid, qcg.to_string()]
#
#                # check to make sure alns match cigar length
#                if aln[4] != qcg.get_len(ref=True):
#                    logger.error(f"Aln length {aln[4]} not matching cigar length on ref {qcg.get_len(ref=True)}! Occurred in {aln}")
#                if aln[5] != qcg.get_len(ref=False):
#                    logger.error(f"Aln length {aln[5]} not matching cigar length on qry {qcg.get_len(ref=False)}! Occurred in {aln}")
#                alns.append(aln)
#
#
#    alns = pd.DataFrame(alns)
#    if alns.empty:
#        return None
#
#    alns = alns.loc[alns[6] > _MIN_SYN_ID] # filter for specified aln identity
#    # count inverted alns as well
#    #alns.loc[alns[8] == -1, 2] = alns.loc[alns[8] == -1, 2] + alns.loc[alns[8] == -1, 3]
#    #alns.loc[alns[8] == -1, 3] = alns.loc[alns[8] == -1, 2] - alns.loc[alns[8] == -1, 3]
#    #alns.loc[alns[8] == -1, 2] = alns.loc[alns[8] == -1, 2] - alns.loc[alns[8] == -1, 3]
#    alns.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr", 'cigar']
#    alns.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
#    return None if alns.empty else alns
#
#
#
#cpdef generate_seqdict(seqh, mappingtrees, chrdict):
#    return {org:('N'*_SPACER_LEN).join([
#        seqh.get_range(Range(org, chrdict[org], interval.data, interval.data + (interval.end - interval.begin)), margin=0)
#        for interval in sorted(mappingtrees[org])])
#        for org in mappingtrees}
#
#


cpdef subtract_mts(mappingtrees, merasyns, skip_ref=True):
    """
    Takes a dict containing an `Intervaltree` with offsets for each organism to be realigned (as produced by `construct_mts`), and returns new mappingtrees with regions covered by multisyn objects in `merasyns` subtracted.
    By default, the reference organism of the merasyn is skipped, as it would be deleted anyway.
    This can be changed by setting `skip_ref` to `False`.
    Used to remove merasynteny found during realignment from the mappingtrees, to do further realignment.
    The initial implementation was just reconstructing the trees at every step, but subtracting is more efficient.
    """
    # core iteration is the same as in construct_mts, except we don't need to store the offsets
    cdef:
        curdict = {org:list(mappingtrees[org][0])[0] for org in mappingtrees if len(mappingtrees[org]) > 0} # stores the current interval in each org
        listdict = defaultdict(list) # used to construct the output mappingtrees
    # these need to be reconstructed to take care of handling the separator intervals

    #print(merasyns)
    #print(type(merasyns))
    for merasyn in merasyns:
        # only subtract on the ref if explicitly specified;
        # would get deleted anyway unless annotating private regions
        for org, rng in\
                merasyn.ranges_dict.items() if skip_ref\
                else [(merasyn.ref.org, merasyn.ref)] + list(merasyn.ranges_dict.items()):
            if not org in curdict or curdict[org] is None: # skip if there are no more intervals to process for this org
                continue
            curint = curdict[org]
            orglist = listdict[org]

            # skip to first interval overlapping this merasyn
            # there should be one of these for any alignment
            while curint.end - curint.begin + curint.data < rng.start:
                orglist.append( (curint.data, curint.end - curint.begin) )
                if curint.end + _SPACER_LEN + 1 in mappingtrees[org]:
                    curint = list(mappingtrees[org][curint.end + _SPACER_LEN + 1])[0]
                else: # there is no offset after this
                    break
            else: # to skip to next org in big for loop
                #logger.warning("Skipped out of int iteration loop!")
                curdict[org] = None
                continue

            # there shouldn't ever be a merasyn spanning beyond one offset
            # directly adjacent offsets are compressed, and other alns should be split
            # emit a warning if this is still the case
            if rng.end > curint.end - curint.begin + curint.data:
                logger.debug(f"{rng.end}, {curint.end - curint.begin + curint.data}")
                logger.warning(f"Synteny in a spacer detected! An alignment went into the separator. Most likely, something went wrong during the alignment call ({rng.end} vs {curint.end - curint.begin + curint.data}).")
            

            # there was no interval overlapping this merasyn anyway, we don't need to subtract anything
            if curint.data > rng.end:
                curdict[org] = curint
                continue
            
            ## from here, rng is fully within curint

            # remove overlap from start, add as separate offset if large enough
            l = rng.start - curint.data
            if l > _MIN_REALIGN_LEN:
                orglist.append( (curint.data, l) )

            # set curint to what remains after removing this merasyn if large enough
            l = curint.data + curint.end - curint.begin - rng.end
            if l > _MIN_REALIGN_LEN:
                curdict[org] = Interval(curint.end - l, curint.end, rng.end)
                #orglist.append( (rng.end, l) )
            else:
                # otherwise skip to next interval
                if curint.end + _SPACER_LEN + 1 in mappingtrees[org]:
                    curdict[org] = list(mappingtrees[org][curint.end + _SPACER_LEN + 1])[0]
                else:
                    curdict[org] = None
                    continue

    return listdict_to_mts(listdict)


cpdef listdict_to_mts(lists):
    """
    Transforms a list of alignments for each organism into a dictionary of intervaltrees mapping an offset in the alignment reference to each position in the sample genome.
    """
    ret = dict()
    posdict = defaultdict(int) # stores the current position in each org
    for org, lst in lists.items():
        tree = IntervalTree()
        for offset, length in lst:
            if length > _MIN_REALIGN_LEN: # filter again for sufficient len
                tree[posdict[org]:posdict[org] + length +1] = offset # end is non-inclusive
                posdict[org] += length + _SPACER_LEN # add interval + spacer

        #if len(tree) > 0:
        ret[org] = tree
            
    return ret

