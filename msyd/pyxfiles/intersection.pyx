#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
#import numpy as np
import functools
from collections import deque
import multiprocessing

#from cython.parallel import prange

import msyd.io
import msyd.util as util
from msyd.orgs import OrgContainer
import msyd.syri_handler as syri_handler
import msyd.cigar
from msyd.multisyn import Multisyn, MultisynContainer, ChromContainer

cdef const int MIN_SYN_THRESH = 50
cdef int SPLIT_INDEL_THRESH = MIN_SYN_THRESH # const as well, but compiler complains when annotating it
cdef volatile int DROPPED_BASES = -1 # global counter, for logging; this isn't perfectly thread safe but that should be fine as we're just logging
# setting to -1 disables logging

logger = util.CustomFormatter.getlogger(__name__)

cpdef int get_SPLIT_INDEL_THRESH():
    global SPLIT_INDEL_THRESH
    return SPLIT_INDEL_THRESH

cpdef set_SPLIT_INDEL_THRESH(int val):
    global SPLIT_INDEL_THRESH
    SPLIT_INDEL_THRESH = val

cpdef int get_dropped_bases():
    global DROPPED_BASES
    return DROPPED_BASES

cpdef start_log_dropped_bases():
    global DROPPED_BASES
    DROPPED_BASES = 0

cpdef stop_log_dropped_bases():
    global DROPPED_BASES
    DROPPED_BASES = -1

cpdef int get_min_syn_thresh():
    return MIN_SYN_THRESH

cdef filter_multisyn(multisyn, drop_small=True, allow_private=False):
    """
    Tests if a multisyn should be added to the output by checking degree and length.
    Unless `drop_private` is set to false, will drop private regions (i.e. merasynteny of degree one)
    If `drop_small` is not set to `False`, will mutate the input multisyn to remove any organisms where the region is smaller than `MIN_SYN_THRESH`.
    """
    global DROPPED_BASES
    if not multisyn: # filter empty objects to handle failures
        return False
    if len(multisyn.ref) < MIN_SYN_THRESH: # filter small regions
        if DROPPED_BASES >= 0: # log dropping if enabled
            DROPPED_BASES += len(multisyn.ref) * multisyn.get_degree()
        return False
    if not multisyn.check(allow_private=allow_private):
        logger.warning(f"{multisyn}")
        return False

    # delete small syntenic regions from the multisyn object, mutates multisyn but that should be fine in this case
    if drop_small:
        droplist = [org for org, rng in multisyn.ranges_dict.items() if len(rng) < MIN_SYN_THRESH]
        for org in droplist:
            if DROPPED_BASES >= 0: # log dropping if enabled
                DROPPED_BASES += len(multisyn.ranges_dict[org])
            del multisyn.ranges_dict[org]
            if multisyn.cigars_dict:
                del multisyn.cigars_dict[org]

    return True
    
cpdef find_overlaps(lmsyncont, rmsyncont, only_core=False, trim=False, allow_private=False):
    """
    This function takes two dataframes containing syntenic regions and outputs the overlap found between each of them as a new pandas dataframe.
    It runs in O(len(left) + len(right)).
    """
    # for finding fully private regions, only_core has to be passed
    #if allow_private:
    #    only_core = True

    cdef:
        object lit = iter(lmsyncont)
        object rit = iter(rmsyncont)
        object ret = MultisynContainer(lmsyncont.orgs + rmsyncont.orgs) # MultisynContainer
        object r # Multisyn, but cython doesn't like this as type annotation
        object l
        int cov = 0 # store the last position in the ref that has been covered in ret
    try:
        r = next(rit)
        l = next(lit)
        cov = 0 # store the last position in the ref that has been covered in ret
    except StopIteration:
        logger.error("Empty iterator passed to find_overlaps!")
        raise ValueError("Empty iterator passed to find_overlaps!")

    while True:
        try: # python iterators suck, so this loop is entirely try-catch'ed
            ovstart = max(r.ref.start, l.ref.start)
            ovend = min(r.ref.end, l.ref.end)

            # find which segment is the starting/ending one
            starting = l if l.ref.start < r.ref.start else r

            if ovend - ovstart >= MIN_SYN_THRESH: # there is valid overlap
                # add the region up to the overlap if it is large enough
                intstart = max(starting.ref.start, cov + 1) # start of the non-overlapping region of interest
                
                if not only_core and ovstart - intstart >= MIN_SYN_THRESH:
                    multisyn = starting.drop(intstart - starting.ref.start, 1 + starting.ref.end - ovstart)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

                # add overlap region
                multisyn = l.drop(ovstart - l.ref.start, l.ref.end - ovend) + r.drop(ovstart - r.ref.start, r.ref.end - ovend)
                if filter_multisyn(multisyn, drop_small=True, allow_private=allow_private):
                    if trim:
                        multisyn.trim_matching_inplace()
                    ret.append(multisyn)
                # everything up to the end of the overlap is covered now
                cov = ovend

            # ratchet by dropping the segment ending first
            # when dropping, include all of the segment that has not been covered by ret so far
            # includes the segment right of an overlap
            if l.ref.end > r.ref.end: # left is after right
                if not only_core and r.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = r.drop(max(0, 1 + cov - r.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

                cov = r.ref.end
                r = next(rit)

            elif r.ref.end > l.ref.end: # right is after left
                if not only_core and l.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = l.drop(max(0, 1 + cov - l.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)
                cov = l.ref.end
                l = next(lit)

            # if they stop at the same position, drop the one starting further left
            elif l.ref.start > r.ref.start:
                if not only_core and r.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = r.drop(max(0, 1 + cov - r.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

                cov = r.ref.end
                r = next(rit)

            else: # do whatever
                if not only_core and l.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = l.drop(max(0, 1 + cov - l.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

                cov = l.ref.end
                l = next(lit)

        except StopIteration: # nothing more to match
            if not only_core and l.ref.chrom == r.ref.chrom: # the loop ended after an overlap call
                # everything up to cov is covered, and starting is guaranteed to be fully covered
                ending = l if l.ref.end > r.ref.end else r
                #print(cov, ending.ref, ending.ref.end - ending.ref.start, {org:cg.get_len(ref=True) for org, cg in ending.cigars_dict.items()})
                if ending.ref.end - cov >= MIN_SYN_THRESH: # check if there is still something to add; if so, add it
                    multisyn = ending.drop(max(0, 1 + cov - ending.ref.start), 0)
                    if filter_multisyn(multisyn, drop_small=True, allow_private=allow_private):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

            break

    if not only_core: # if calling crosssyn, also add remaining multisyn if there is any
        for l in lit:
            if filter_multisyn(l):
                if trim:
                    l.trim_matching_inplace()
                ret.append(l)
        for r in rit:
            if filter_multisyn(r):
                if trim:
                    r.trim_matching_inplace()
                ret.append(r)

    if len(ret) == 0:
        logger.error("find_multisynteny found no overlapping synteny!")
    return ret # shouldn't need sorting
#END

cpdef object split_indels(msyncont: MultisynContainer):
    cdef object ret = MultisynContainer(orgs = msyncont.orgs, cap=len(msyncont))
    # split the multisyns if there are any large indels in the alignments
    for multisyn in iter(msyncont):
        # filter out short multisyns here
        ret.extend([msyn for msyn in multisyn.split_indels(SPLIT_INDEL_THRESH)\
                if filter_multisyn(msyn)]) # [] to work around cython not liking lambdas
    return ret

cpdef find_multisyn(qrynames, syris, alns, cores=1, base=None, sort=False, ref='a', refname="ref", SYNAL=True, disable_overlapcheck=False, only_core=False, trim=True):
    """
    Finds core and cross-syntenic regions containing the reference in the input files, depending on if the parameter `only_core` is `True` or `False`.
    Fairly conservative.
    Uses either SYNAL or SYN regions as annotated by SyRI, controlled by the parameter `SYNAL`.
    In the case of SYN regions, alignment-based length calculation is not supported and `alns` is ignored.

    :param: a list of filenames of SyRI output and alignment files in BAM, SAM or PAF format to read in, parameters
    optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be
     sorted (default False).
    `alns` can be set to `None`, in which case the region lengths will be estimated instead of calculated exactly from CIGAR strings.
    :param only_core: Whether to output all cross synteny or only core syntenic regions.
    :return: a pandas dataframe containing the chromosome, start and end positions of the core syntenic region for each organism.
    """

    chromsyn = prepare_input(qrynames, syris, alns, cores=cores, base=base, sort=sort, ref=ref, SYNAL=SYNAL)
    logger.info("Finished reading input files, starting intersection.")

    #NOTE pass along cores/len(chroms) as cores to here?
    _process_synlists = functools.partial(process_synlists, refname=refname, only_core=only_core, trim=trim)
    return chromsyn.apply_chroms_par(_process_synlists, ncores=cores)


cpdef prepare_input(qrynames, syris, alns, refname="ref", cores=1, base=None, sort=False, ref='a', SYNAL=True, disable_overlapcheck=False):
    """
    Fetches input from filenames given to it; mostly parallelized.
    :Returns: a Dict of chromosome IDs to a list of Multisyn DFs (one per sample).
    This allows seamless parallelization between chromosome IDs
    """
    from msyd.syri_handler import extract_from_filelist, match_synal

    ## read in as dict from chrnames to a list containing a dataframe with all syri SYNAL calls for each input org file
    orgs = OrgContainer.from_list(qrynames)
    orgs.add_org(refname) # add ref to orgs
    syndict = extract_from_filelist(syris, qrynames, cores=cores, anns=["SYNAL"] if SYNAL else ["SYN"])
    if sort:
        syndict = {chrom: [syndf.sort_values(syndf.columns[0]) for syndf in syndfs]for chrom, syndfs in syndict}

    #alnfilelookup = {
    #        'sam': msyd.io.readSAMBAM,
    #        'bam': msyd.io.readSAMBAM,
    #        'paf': msyd.io.readPAF
    #        }

    ## non-cigar path
    ## this code path shouldn't really be used anymore.
    if not (SYNAL and alns):
        logger.warning("No alignments found or `--syn` passed! Assuming all synteny to be exactly identical. This is fast but error-prone and inaccurate.")
        return ChromContainer({chrom:[pd.DataFrame([Multisyn(ref=row[1][0],
                        ranges_dict={row[1][1].org:row[1][1]}, cigars_dict = None)
                        for row in s.iterrows()]) for s in syns]
                for chrom, syns in syndict}, orgs)

    #with multiprocessing.Pool(cores) as pool:
    #    alns = pool.map(lambda aln: msyd.io.alnfilelookup[aln.split('.')[-1]](aln), alns)
    #    alns = pool.map(lambda aln: aln[(aln.adir==1) & (aln.bdir==1)], alns) # pre-filter to non-inverted alns
    alns = [msyd.io.alnfilelookup[aln.split('.')[-1]](aln) for aln in alns]
    alns = [aln[(aln.adir==1) & (aln.bdir==1)] for aln in alns] # pre-filter to non-inverted alns
        
    # this step is single-threaded; TODO parallelize?
    alndict = msyd.io.collate_by_chrom(alns, chromid=ref+'chr')

    for chrom in syndict: #TODO maybe parallelize over chrs instead
        #syndict[chrom] = pool.map(lambda syndf, alndf: match_synal(syndf, alndf, ref=ref), zip(syndict[chrom], alndict[chrom]))
        syndict[chrom] = [match_synal(syndf, alndf, ref=ref, refname=refname) for syndf, alndf in zip(syndict[chrom], alndict[chrom])]
    
    return ChromContainer(dict(syndict), orgs)

cpdef process_synlists(chrom, msyncontlist, base=None, disable_overlapcheck=False, cores=1, only_core=False, trim=True):
    # msyncontlist is a list of MultisynContainers

    # remove overlap
    #NOTE could multithread this, but not sure its worth it
    if not disable_overlapcheck:
        for msyncont in msyncontlist:
            syri_handler.handle_conflicts(iter(msyncont))

    if SPLIT_INDEL_THRESH > 0:
        logger.info(f"{chrom}: Splitting alignments at indels > {SPLIT_INDEL_THRESH} bp")
        if cores == 1:
            msyncontlist = [split_indels(msyncont) for msyncont in msyncontlist]
        else:
            with multiprocessing.Pool(cores) as pool:
                msyncontlist = pool.map(split_indels, msyncontlist)

    logger.info(f"{chrom}: overlapping synteny trimmed")

    # shouldn't need any overlap removal
    if base:
        logger.info("Adding msyns from PSF file for incremental calling")
        msyncontlist.append(base)

    return reduce_find_overlaps(msyncontlist, cores, only_core=only_core, trim=trim)
# END

cpdef reduce_find_overlaps(syns, cores, only_core=False, trim=True):
    if len(syns) == 0:
        return None
    multisyns = None
    ovlap = functools.partial(find_overlaps, only_core=only_core, trim=trim)
    if cores > 1:
        multisyns = util.parallel_reduce(ovlap, syns, cores)
    else:
        multisyns = functools.reduce(ovlap, syns)

    return multisyns


