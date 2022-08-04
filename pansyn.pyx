#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
#import numpy as np
import copy
import functools
from collections import deque

import ingest
import util
from cigar import Cigar
from classes import Pansyn, Range, Position

MIN_SYN_THRESH = 0


def calc_overlap(l: Pansyn, r: Pansyn, detect_crosssyn=False, allow_overlap=False):
    """
    Used to be called `combine` in coresyn and crosssyn.
    Takes two `Pansyn` objects, finds the overlap and returns a list of the `Pansyn` objects generated by overlapping the two.
    Returns a sorted list of at most three `Pansyn` objects.
    :param: `detect_crosssyn` specifies whether to output anything but the overlap.
    If disabled will only find core syntenic regions.
    `allow_overlap` specifies whether pansyntenic regions are allowed to overlap, in essence outputting the input regions as-is.
    Output may be non-associative and non-commutative if set to `True`.
    """

    if not l.ref and not r.ref:
        raise ValueError("ERROR: calc_overlap can only be called on two Pansyns with references!")

    # compute overlap
    ovstart = max(l.ref.start, r.ref.start)
    ovend = min(l.ref.end, r.ref.end)
    if not ovstart < ovend:
        print(f"WARNING: no overlap found between {l} and {r}, returning empty!")
        return []

    leftest = l if l.ref.start < r.ref.start else r
    rightest = l if l.ref.end > r.ref.end else r

    ret = set() # use a set to automatically remove duplicates

    def add_filtered(pansyn): # checks if the region is higher than MIN_SYN_THRESH, if so, adds it to ret
        if not pansyn:
            return
        if len(pansyn.ref) < MIN_SYN_THRESH:
            return
        if pansyn.get_degree() < 1:
            return
        if not detect_crosssyn:
            if pansyn.get_degree() < l.get_degree() + r.get_degree():
                return
        ret.add(pansyn)



    if detect_crosssyn:
        if allow_overlap:
            add_filtered(leftest)
        else:
            add_filtered(leftest.drop(0, leftest.ref.end - ovstart))
    
    # core synteny
    add_filtered(l.drop(ovstart - l.ref.start, l.ref.end - ovend) + r.drop(ovstart - r.ref.start, r.ref.end - ovend))

    if detect_crosssyn:
        if allow_overlap:
            add_filtered(rightest)
        else:
            add_filtered(rightest.drop(rightest.ref.start - ovend, 0))

    return sorted(ret)


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
                cg = Cigar.from_string(alnr['cg'])
                rng = synr[1]
                pansyn = Pansyn(ref=synr[0], ranges_dict={org:rng}, cigars_dict={org:cg})
                #rng.end = rng.start + cg.get_len(ref=False) -1
                ## debugging output
                #try:
                #    pansyn.check()
                #except ValueError as  e:
                #    print(e)
                #    print(f"ERRORERRROR: cigar string length {cg.get_len(ref=True)}, {cg.get_len(ref=False)} does not match ref or qry length {len(pansyn.ref)}/{len(list(pansyn.ranges_dict.values())[0])} in {pansyn}")
                #    print(cg)
                #    print(pansyn)
                #    print("Adjusting to cg length")
                #    # forcibly ajust end position to match cigar length, as that doesn't always seem to be the case in syri/pysam output for some reason
                #    rng.end = rng.start + cg.get_len(ref=False) -1

                ret.append(pansyn)
                synr = next(syniter)[1]
            alnr = next(alniter)[1]
        except StopIteration:
            break

    return pd.DataFrame(list(ret))


def remove_overlap(syn):
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
        # check for overlap on the reference
        ov = prev.ref.end - cur.ref.start
        if ov <= 0 or cur.ref.chr != prev.ref.chr: # there is no overlap
            prev = cur
            continue

        #print(ov, prevref, curref)
        # remove the overlap from the previous row;
        # this performs essentially the same function as compute_overlap
        # in intersect_syns
        cur.ref.start += ov
        if cur.cigars_dict:
            for org, cg in cur.cigars_dict.items():
                drop, cgnew = cg.get_removed(ov, start=True)
                cur.cigars_dict[org] = cgnew
                cur.ranges_dict[org].start += drop
        else:
            for rng in cur.ranges_dict.values():
                rng.start += ov

        prev = cur


def find_overlaps(left, right, **kwargs):
    """
    This function takes two dataframes containing syntenic regions and outputs the overlap found between each of them as a new pandas dataframe.
    It runs in O(len(left) + len(right)).
    Calls `calc_overlap` to determine the actual overlaps to merge as output.
    """
    ret = deque()

    if len(right) == 0:
        raise ValueError("right is empty!")
    if len(left) == 0:
        raise ValueError("left is empty!")

    riter = right.iterrows()
    rrow = next(riter)[1][0]
    liter = left.iterrows()
    lrow = next(liter)[1][0]
    while True:
        try: # python iterators suck, so this loop is entirely try-catch'ed
            lref = lrow.ref
            rref = rrow.ref
            if rref.chr > lref.chr:
                lrow = next(liter)[1][0]
                continue
            if lref.chr > rref.chr:
                rrow = next(riter)[1][0]
                continue
            
            # determine if there is an overlap
            ovstart = max(rref.start, lref.start)
            ovend = min(rref.end, lref.end)
            if ovend - ovstart > MIN_SYN_THRESH: # there is valid overlap
                ret.extend(calc_overlap(lrow, rrow, **kwargs))

            # ratchet by dropping the segment with a smaller end
            if lref.end > rref.end: # left is after right
                rrow = next(riter)[1][0]
            elif rref.end > lref.end: # right is after left
                lrow = next(liter)[1][0]
                # if they stop at the same position, drop the one starting further left
            elif lref.start > rref.start:
                rrow = next(riter)[1][0]
            else: # do whatever
                lrow = next(liter)[1][0]

        except StopIteration: # nothing more to match
            break

    del riter
    del liter
    # sort to be safe, maybe optimitze this away later?
    ret = pd.DataFrame(data=list(ret))
    return ret.sort_values(ret.columns[0])


def find_multisyn(syris, alns, sort=False, ref='a', cores=1, SYNAL=True, **kwargs):
    """
    Finds core and cross-syntenic regions in the input files, depending on if the parameter `detect_crossyn` that is ultimately passed on to `calc_overlap` is set to `True`.
    Fairly conservative.
    :param: a list of filenames of SyRI output and alignment files in BAM, SAM or PAF format to read in, parameters optionally specifying which sequence is the reference (default 'a') and a boolean specifying if the input needs to be sorted (default False).
    `alns` can be set to `None`, in which case the region lengts will be estimated instead of calculated exactly from CIGAR strings.
    Crosssynteny is detected if the `detect_crosssyn` parameter is set to `True`.
    :return: a pandas dataframe containing the chromosome, start and end positions of the core syntenic region for each organism.
    """

    syns = util.extract_regions_to_list(syris, ann="SYNAL" if SYNAL else "SYN")

    if sort:
        syns = [x.sort_values(x.columns[0]) for x in syns]

    alnfilelookup = {
            'sam': ingest.readSAMBAM,
            'bam': ingest.readSAMBAM,
            'paf': ingest.readPAF
            }

    if alns and SYNAL:
        #TODO log
        alns = [alnfilelookup[aln.split('.')[-1]](aln) for aln in alns]
        alns = [aln[(aln.adir==1) & (aln.bdir==1)] for aln in alns] # only count non-inverted alignments as syntenic

        syns = [match_synal(*x, ref=ref) for x in zip(syns, alns)]
    else:
        syns = [
                pd.DataFrame(
                    [ Pansyn(ref=row[1][0],
                        ranges_dict={row[1][1].org:row[1][1]}, cigars_dict = None)
                        for row in s.iterrows()]) for s in syns]

    # remove overlap
    for syn in syns:
        remove_overlap(syn)
    print("INFO: overlap removed")

    pansyns = None
    ovlap = functools.partial(find_overlaps, **kwargs)
    if cores > 1:
        pansyns = util.parallel_reduce(ovlap, syns, cores)
    else:
        pansyns = functools.reduce(ovlap, syns)

    return pansyns

