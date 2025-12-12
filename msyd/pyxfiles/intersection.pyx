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

import msyd.io as io
import msyd.util as util
import msyd.syri_handler as syri_handler
import msyd.cigar
from msyd.multisyn import Multisyn

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
    
cpdef find_overlaps(left, right, only_core=False, trim=False, allow_private=False):
    """
    This function takes two dataframes containing syntenic regions and outputs the overlap found between each of them as a new pandas dataframe.
    It runs in O(len(left) + len(right)).
    """
    cdef ret = deque()
    #print("l:", left, "r:", right)

    # for finding fully private regions, only_core has to be passed
    #if allow_private:
    #    only_core = True

    if right is None or right.empty:
        logger.error("find_overlap called with no multisyn regions (right)!")
        #raise ValueError("right is empty!")
        return None
    if left is None or left.empty:
        logger.error("find_overlap called with no multisyn regions (left)!")
        #raise ValueError("left is empty!")
        return None

    cdef:
        rit = right.iterrows()
        lit = left.iterrows()
        r = next(rit)[1][0]
        l = next(lit)[1][0]

    cdef int cov = 0 # store the last position in the ref that has been covered in ret

    ## helper Fn to only add filtered multisyns to the final output
    ## calls filter_multisyn and if only_core is set additionally filters for core synteny
    #def add_filtered(multisyn):
    #    if filter_multisyn(multisyn):
    #        # filter non-core-syntenic regions in case that's relevant
    #        # don't double count reference
    #        if only_core and multisyn.get_degree() < (l.get_degree() + r.get_degree() - 1):
    #            return
    #        ret.append(multisyn)
    
    while True:
        try: # python iterators suck, so this loop is entirely try-catch'ed
            # ensure that the chr matches, reset the covered region
            if r.ref.chr > l.ref.chr:
                cov = -1
                l = next(lit)[1][0]
                continue
            if l.ref.chr > r.ref.chr:
                cov = -1
                r = next(rit)[1][0]
                continue
            
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
                r = next(rit)[1][0]

            elif r.ref.end > l.ref.end: # right is after left
                if not only_core and l.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = l.drop(max(0, 1 + cov - l.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)
                cov = l.ref.end
                l = next(lit)[1][0]

            # if they stop at the same position, drop the one starting further left
            elif l.ref.start > r.ref.start:
                if not only_core and r.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = r.drop(max(0, 1 + cov - r.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

                cov = r.ref.end
                r = next(rit)[1][0]

            else: # do whatever
                if not only_core and l.ref.end - cov >= MIN_SYN_THRESH:
                    multisyn = l.drop(max(0, 1 + cov - l.ref.start), 0)
                    if filter_multisyn(multisyn):
                        if trim:
                            multisyn.trim_matching_inplace()
                        ret.append(multisyn)

                cov = l.ref.end
                l = next(lit)[1][0]

        except StopIteration: # nothing more to match
            if not only_core and l.ref.chr == r.ref.chr: # the loop ended after an overlap call
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
            l = l[1][0]
            if filter_multisyn(l):
                if trim:
                    l.trim_matching_inplace()
                ret.append(l)
        for r in rit:
            r = r[1][0]
            if filter_multisyn(r):
                if trim:
                    r.trim_matching_inplace()
                ret.append(r)

    if len(ret) == 0:
        return pd.DataFrame()
    else:
        return pd.DataFrame(data=list(ret)) # shouldn't need sorting
#END

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


cpdef split_indels(syndf):
    # split the multisyns if there are any large indels in the alignments
    cdef ret = deque()

    for _, multisyn in syndf.iterrows():
        multisyn = multisyn[0]
        # filter out short multisyns here
        ret.extend([msyn for msyn in multisyn.split_indels(SPLIT_INDEL_THRESH)\
                if filter_multisyn(msyn)])
        #for msyn in multisyn.split_indels(SPLIT_INDEL_THRESH):
        #    if filter_multisyn(msyn):
        #        ret.append(msyn)
        #ret.extend(filter(filter_multisyn, multisyn.split_indels(SPLIT_INDEL_THRESH))) # doesn't work for some reason?
    return pd.DataFrame(list(ret))


cpdef find_multisyn(qrynames, syris, alns, cores=1, base=None, sort=False, ref='a', SYNAL=True, disable_overlapcheck=False, only_core=False, trim=True):
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

    syndict = prepare_input(qrynames, syris, alns, cores=cores, base=base, sort=sort, ref=ref, SYNAL=SYNAL)
    logger.info("Finished reading input files, starting intersection.")

    return process_syndicts(syndict, cores=cores, only_core=only_core, trim=trim)


def process_syndicts(syndict, cores=4, only_core=False, trim=True):
    """
    Small fn to do parallel processing of a dictionary of syndfs per chromosome.
    """
    ## NOTE THOUGHTS:
    ## make syndict into ChromContainer
    ## use mapreduce call to iter_chrom_par to get rid of this fn
    ## should simplify some code
    if not syndict:
        return {}

    _process_syndf = functools.partial(process_syndfs, only_core=only_core)
    chroms, syndfs = zip(*syndict.items())

    if cores > 1:
        with multiprocessing.Pool(cores) as pool:
            results = pool.map(_process_syndf, syndfs)
    else:
        results = map(_process_syndf, syndfs)
    # map guarantees the order of results is the same as the order of input
    return dict(zip(chroms, results))

    #cdef list chromlist = list(syndict)
    #cdef int n = len(chromlist)
    #cdef int i
    #for i in prange(n, nogil=True):
    #    with gil:
    #        chrom = chromlist[i]
    #        syndf = syndict[chrom]
    #    intersected = process_syndfs(syndf)
    #    with gil:
    #        syndict[chrom] = intersected
    # return syndict

cpdef _workaround(tup): # tup: [chrom, syndfs, dict kwargs]
    # Annoying workaround, because multiprocessing doesn't like lambdas
    return tup[0], process_syndfs(tup[1], **tup[2])


cpdef prepare_input(qrynames, syris, alns, cores=1, base=None, sort=False, ref='a', SYNAL=True, disable_overlapcheck=False):
    """
    Fetches input from filenames given to it; mostly parallelized.
    :Returns: a Dict of chromosome IDs to a list of Multisyn DFs (one per sample).
    This allows seamless parallelization between chromosome IDs
    """
    from msyd.syri_handler import extract_from_filelist, match_synal

    syndict = extract_from_filelist(syris, qrynames, cores=cores, anns=["SYNAL"] if SYNAL else ["SYN"])
    if sort:
        syndict = {chrom: [syndf.sort_values(syndf.columns[0]) for syndf in syndfs]for chrom, syndfs in syndict}

    #alnfilelookup = {
    #        'sam': io.readSAMBAM,
    #        'bam': io.readSAMBAM,
    #        'paf': io.readPAF
    #        }

    if not (SYNAL and alns):
        logger.warning("No alignments found or `--syn` passed! Assuming all synteny to be exactly identical. This is fast but error-prone and inaccurate.")
        return {chrom:[pd.DataFrame([Multisyn(ref=row[1][0],
                        ranges_dict={row[1][1].org:row[1][1]}, cigars_dict = None)
                        for row in s.iterrows()]) for s in syns]
                for chrom, syns in syndict}

    #with multiprocessing.Pool(cores) as pool:
    #    alns = pool.map(lambda aln: io.alnfilelookup[aln.split('.')[-1]](aln), alns)
    #    alns = pool.map(lambda aln: aln[(aln.adir==1) & (aln.bdir==1)], alns) # pre-filter to non-inverted alns
    alns = [io.alnfilelookup[aln.split('.')[-1]](aln) for aln in alns]
    alns = [aln[(aln.adir==1) & (aln.bdir==1)] for aln in alns] # pre-filter to non-inverted alns
        
    # this step is single-threaded; TODO parallelize?
    alndict = io.collate_by_chrom(alns, chromid=ref+'chr')

    for chrom in syndict: #TODO maybe parallelize over chrs instead
        #syndict[chrom] = pool.map(lambda syndf, alndf: match_synal(syndf, alndf, ref=ref), zip(syndict[chrom], alndict[chrom]))
        syndict[chrom] = [match_synal(syndf, alndf, ref=ref) for syndf, alndf in zip(syndict[chrom], alndict[chrom])]
            
    return syndict

cpdef process_syndfs(syndfs, base=None, disable_overlapcheck=False, cores=1, only_core=False, trim=True):
    # remove overlap
    if not disable_overlapcheck:
        if cores == 1:
            syndfs = [syri_handler.handle_conflicts(syndf) for syndf in syndfs]
        else:
            with multiprocessing.Pool(cores) as pool:
                syndfs = pool.map(syri_handler.handle_conflicts, syndfs)

    if SPLIT_INDEL_THRESH > 0:
        logger.info(f"Splitting alignments at indels > {SPLIT_INDEL_THRESH} bp")
        if cores == 1:
            syndfs = [split_indels(syndf) for syndf in syndfs]
        else:
            with multiprocessing.Pool(cores) as pool:
                syndfs = pool.map(split_indels, syndfs)

    logger.info("overlapping synteny trimmed")

    # shouldn't need any overlap removal
    if base:
        logger.info("reading in PSF for incremental calling")
        syndfs.append(base)

    return reduce_find_overlaps(syndfs, cores, only_core=only_core, trim=trim)
# END
