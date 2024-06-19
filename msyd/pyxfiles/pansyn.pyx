#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
#import numpy as np
import copy
import functools
from collections import deque
import multiprocessing

import msyd.io as io
import msyd.scripts.util as util
import msyd.cigar
from msyd.cigar import Cigar
from msyd.coords import Pansyn, Range, Position

cdef int MIN_SYN_THRESH = 30

logger = util.CustomFormatter.getlogger(__name__)
import logging
logger.setLevel(logging.DEBUG)

def find_overlaps(left, right, only_core=False):
    """
    This function takes two dataframes containing syntenic regions and outputs the overlap found between each of them as a new pandas dataframe.
    It runs in O(len(left) + len(right)).
    """
    ret = deque()
    #print("l:", left, "r:", right)

    if right is None or right.empty:
        logger.error("find_overlap called with no pansyn regions (right)!")
        #raise ValueError("right is empty!")
        return None
    if left is None or left.empty:
        logger.error("find_overlap called with no pansyn regions (left)!")
        #raise ValueError("left is empty!")
        return None

    rit = right.iterrows()
    lit = left.iterrows()
    r = next(rit)[1][0]
    l = next(lit)[1][0]

    cdef int cov = 0 # store the last position in the ref that has been covered in ret

    #TODO refactor into external fn, write filtering fn for reduce_find_overlap if n=1
    def add_filtered(pansyn): # Filters out pansyns that should be discarded at this step, and adds those that we want to keep to the results list
        if not pansyn: # filter empty objects to handle failures
            return
        if len(pansyn.ref) < MIN_SYN_THRESH: # filter small regions
            return

        # delete small syntenic regions from the pansyn object, mutates pansyn but that should be fine in this case
        remlist = [org for org, rng in pansyn.ranges_dict.items() if len(rng) < MIN_SYN_THRESH]
        for org in remlist:
            del pansyn.ranges_dict[org]
            if pansyn.cigars_dict:
                del pansyn.cigars_dict[org]

        if pansyn.get_degree() < 1: # filter regions without non-reference synteny
            return
        if only_core:
            if pansyn.get_degree() < l.get_degree() + r.get_degree(): # filter non-core-syntenic regions in case that's relevant
                return
        ret.append(pansyn)
    
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
                    add_filtered(starting.drop(intstart - starting.ref.start, 1 + starting.ref.end - ovstart))

                # add overlap region
                add_filtered(l.drop(ovstart - l.ref.start, l.ref.end - ovend) + r.drop(ovstart - r.ref.start, r.ref.end - ovend))
                # everything up to the end of the overlap is covered now
                cov = ovend

            # ratchet by dropping the segment ending first
            # when dropping, include all of the segment that has not been covered by ret so far
            # includes the segment right of an overlap
            if l.ref.end > r.ref.end: # left is after right
                if not only_core and r.ref.end - cov >= MIN_SYN_THRESH:
                    add_filtered(r.drop(max(0, 1 + cov - r.ref.start), 0))
                cov = r.ref.end
                r = next(rit)[1][0]

            elif r.ref.end > l.ref.end: # right is after left
                if not only_core and l.ref.end - cov >= MIN_SYN_THRESH:
                    add_filtered(l.drop(max(0, 1 + cov - l.ref.start), 0))
                cov = l.ref.end
                l = next(lit)[1][0]

                # if they stop at the same position, drop the one starting further left
            elif l.ref.start > r.ref.start:
                if not only_core and r.ref.end - cov >= MIN_SYN_THRESH:
                    add_filtered(r.drop(max(0, 1 + cov - r.ref.start), 0))
                cov = r.ref.end
                r = next(rit)[1][0]

            else: # do whatever
                if not only_core and l.ref.end - cov >= MIN_SYN_THRESH:
                    add_filtered(l.drop(max(0, 1 + cov - l.ref.start), 0))
                cov = l.ref.end
                l = next(lit)[1][0]

        except StopIteration: # nothing more to match
            if not only_core and l.ref.chr == r.ref.chr: # the loop ended after an overlap call
                # everything up to cov is covered, and starting is guaranteed to be fully covered
                ending = l if l.ref.end > r.ref.end else r
                #print(cov, ending.ref, ending.ref.end - ending.ref.start, {org:cg.get_len(ref=True) for org, cg in ending.cigars_dict.items()})
                if ending.ref.end - cov >= MIN_SYN_THRESH: # check if there is still something to add; if so, add it
                    add_filtered(ending.drop(max(0, 1 + cov - ending.ref.start), 0))

            break

    if not only_core: # if calling crosssyn, also add remaining pansyn if there is any
        for l in lit:
            l = l[1][0]
            add_filtered(l)
        for r in rit:
            r = r[1][0]
            add_filtered(l)

    del rit
    del lit
    if len(ret) == 0:
        return pd.DataFrame()
    
    ret = pd.DataFrame(data=list(ret))#sorted(list(ret))) # sorting shouldn't be necessary

    total_len_left = sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], left.iterrows())))
    total_len_right = sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], right.iterrows())))
    total_len_ret = sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], ret.iterrows())))
    logger.debug(f"left orgs: {util.get_orgs_from_df(left)}, right orgs: {util.get_orgs_from_df(right)}, ret orgs: {util.get_orgs_from_df(ret)}")
    logger.debug(f"left len: {total_len_left}, right len: {total_len_right}, ret len: {total_len_ret}")

    return ret#.sort_values(ret.columns[0])
#END


# given a bam file and corresponding SYNAL range df,
# Transform them into one list of Pansyn objects
def match_synal(syn, aln, ref='a'):
    """
    This function takes an aligment and SYNAL dataframe and matches corresponding regions.
    It returns a dataframe containing the regions with the corresponding CIGAR string as a `Pansyn` object.
    :params: syn: SYNAL dataframe, aln: alignment dataframe, ref: whether the reference is the 'a' or 'b' strand in the alignment dataframe.
    :returns: a dataframe containing the SYNAL regions with corresponding CIGAR strings as `Pansyn` objects.
    """
    ret = deque()
    syniter = syn.iterrows()
    alniter = aln.iterrows()
    refchr = ref + "chr"
    refstart = ref + "start"
    refend = ref + "end"

    synr = next(syniter)[1]
    alnr = next(alniter)[1]
    while True:
        try:
            org = synr[1].org
            if synr[0].chr == alnr[refchr] and synr[0].start == alnr[refstart] and synr[0].end == alnr[refend]:
                cg = msyd.cigar.cigar_from_string(alnr['cg'])
                rng = synr[1]
                pansyn = Pansyn(ref=synr[0], ranges_dict={org:rng}, cigars_dict={org:cg})
                #rng.end = rng.start + cg.get_len(ref=False) -1

                ## Correct mismatches between CIGAR and coordinate len

                # check on ref
                if not len(pansyn.ref) == cg.get_len():
                    logger.warning(f"CIGAR len ({cg.get_len()}) not matching coordinate len ({len(pansyn.ref)}) on ref! Adjusting end to match CIGAR len (this might be because of clipping).")
                    pansyn.ref.end = pansyn.ref.start + cg.get_len() - 1

                # check on org
                if not len(rng) == cg.get_len(ref=False):
                    # forcibly ajust end position to match cigar length, as that doesn't always seem to be the case in syri/pysam output for some reason
                    logger.warning(f"CIGAR len ({cg.get_len(ref=False)}) not matching coordinate len ({len(rng)}) on {org}! Adjusting end to match CIGAR len (this might be because of clipping).")
                    rng.end = rng.start + cg.get_len(ref=False) -1


                ret.append(pansyn)
                synr = next(syniter)[1]
            alnr = next(alniter)[1]
        except StopIteration:
            break
    return pd.DataFrame(list(ret))


cdef remove_overlap(syn):
    """
    part of the preprocessing of SYNAL regions for find_multisyn
    removes overlap from the first region if two overlapping regions are next to each other
    assumes syn to be sorted
    mutates syn
    """
    syniter = syn.iterrows()
    prev = next(syniter)[1][0]
    for _, cur in syniter:
        cur = cur[0]
        if cur.ref.chr != prev.ref.chr: # there can be no overlap between chrs
            prev = cur
            continue

        ## check for & remove overlap on the reference
        ov = prev.ref.end - cur.ref.start +1
        if ov > 0:
            # there is overlap on ref
            logger.warning(f"Found {ov} bp overlap on reference at {cur.ref.start}, dropping from latter record!")
            cur.drop_inplace(ov, 0) # call drop_inplace to mutate the dataframe from a reference

        ## check for overlap on other orgs

        # when this is called, cur and prev should normally have the same orgs
        # will not catch overlap between non-adjacent regions!
        #assert(set(cur.ranges_dict) == set(prev.ranges_dict))
        for org in cur.ranges_dict: # should be on the same chr
            if org not in prev.ranges_dict:
                continue
            assert(cur.ranges_dict[org].chr == prev.ranges_dict[org].chr)

            ov = prev.ranges_dict[org].end - cur.ranges_dict[org].start + 1 # indices are inclusive
            if ov > 0:
                # check if the region is fully contained, in case this ever happens
                # drop the region on this org in that case
                if cur.ranges_dict[org].end <= prev.ranges_dict[org].end:
                    del cur.ranges_dict[org]
                    if cur.cigars_dict:
                        del cur.cigars_dict[org]

                # there is overlap on org
                logger.warning(f"Found {ov} bp overlap on {org} at {cur.ranges_dict[org].start}, dropping from latter record!")
                logger.debug(f"Overlapping on {org}: {prev}, {cur}")
                #logger.debug(f"Before drop: {cur}")
                cur.drop_on_org_inplace(ov, 0, org)
                #logger.debug(f"After drop: {cur}")

        prev = cur

    return syn
# END

def find_multisyn(qrynames, syris, alns, base=None, sort=False, ref='a', cores=1, SYNAL=True, overlapping=True, **kwargs):
    """
    Finds core and cross-syntenic regions containing the reference in the input files, depending on if the parameter `only_core` is `True` or `False`.
    Fairly conservative.
    Uses either SYNAL or SYN regions as annotated by SyRI, controlled by the parameter `SYNAL`.
    In the case of SYN regions, alignment-based length calculation is not (yet) supported and `alns` is ignored.

    :param: a list of filenames of SyRI output and alignment files in BAM, SAM or PAF format to read in, parameters
    optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be
     sorted (default False).
    `alns` can be set to `None`, in which case the region lengths will be estimated instead of calculated exactly from CIGAR strings.
    :param only_core: Whether to output all cross synteny or only core syntenic regions.
    :return: a pandas dataframe containing the chromosome, start and end positions of the core syntenic region for each organism.
    """
    from msyd.scripts.io import extract_syri_regions_to_list_from_files

    syns = extract_syri_regions_to_list_from_files(syris, qrynames, cores=cores, anns=["SYNAL"] if SYNAL else ["SYN"])
    if sort:
        syns = [x.sort_values(x.columns[0]) for x in syns]

    #alnfilelookup = {
    #        'sam': io.readSAMBAM,
    #        'bam': io.readSAMBAM,
    #        'paf': io.readPAF
    #        }

    if alns and SYNAL:
        #TODO log
        alns = [io.alnfilelookup[aln.split('.')[-1]](aln) for aln in alns]
        alns = [aln[(aln.adir==1) & (aln.bdir==1)] for aln in alns] # only count non-inverted alignments as syntenic

        syns = [match_synal(*x, ref=ref) for x in zip(syns, alns)]
    else:
        syns = [
                pd.DataFrame(
                    [ Pansyn(ref=row[1][0],
                        ranges_dict={row[1][1].org:row[1][1]}, cigars_dict = None)
                        for row in s.iterrows()]) for s in syns]

    # remove overlap
    if overlapping:
        if cores == 1:
            syns = [remove_overlap(syn) for syn in syns]
        else:
            with multiprocessing.Pool(cores) as pool:
                syns = pool.map(remove_overlap, syns)

    logger.info("overlap removed")

    # shouldn't need any overlap removal
    if base:
        logger.info("reading in PFF for incremental calling")
        syns.append(io.read_pff(base))

    return reduce_find_overlaps(syns, cores, **kwargs)
# END

def reduce_find_overlaps(syns, cores, **kwargs):
    if len(syns) == 0:
        return None
    pansyns = None
    ovlap = functools.partial(find_overlaps, **kwargs)
    if cores > 1:
        pansyns = util.parallel_reduce(ovlap, syns, cores)
    else:
        pansyns = functools.reduce(ovlap, syns)

    return pansyns
