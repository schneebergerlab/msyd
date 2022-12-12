#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

from pansyri.classes.cigar import Cigar#, cig_clips, cig_aln_types
import pansyri.classes as classes

# redeclare, as cdefs can apparently not be imported
cdef cig_types = set(['M', '=', 'X', 'S', 'H', 'D', 'I', 'N'])
cdef cig_aln_types = set(['M', 'X', '='])
cdef cig_clips = set(['S', 'H', 'P', 'N']) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these
import logging

logger = logging.getLogger(__name__)


"""
TODOs
– write unit tests
– write frontend:
    • maybe have a convenience function imputing everything given two bam files?
    • maybe have a convenience function automatically detecting core/cross synteny from such imputed bams?
    • maybe have a fucntion that automatically imputes along all pansyntenic regions?
    => convenience functions may be moved to util
"""

def impute_strings(strl: str, strr: str):
    """Convenience function performing the imputation on just two CIGAR strings.

    :param strl, strr: the two strings to impute the alignment from.
    :return: A CIGAR string containing an approximated alignment of `strr` to `strl`
    :rtype: str
    """
    cgl, cgr = classes.cigar.cigar_from_string(strl), Cigar.from_string(strr)
    return impute(cgl, cgr).to_string()


def impute_ranges(cgl, rngl, cgr, rngr):
    """A function for imputing two regions not corresponding 1:1 in start/end position by padding appropriately.
    Emits a warning if the two regions do not have an overlap.

    :param cgl, cgr: Two `Cigar` objects containing the alignments to the common reference for each of the two positions.
    :param rngl, rngr: Two `Range` objects containing the positions on the common reference
    :return: A `Cigar` containing the approximated alignment of the two `Cigars`.
    :rtype: `Cigar`
    """#TODO maybe also compute Range of alignment on target alignment? issue: needs to supply starting point in methot arguments => move to different method or use default arguments to not do this as standard
    maxstart = max(rngl.start, rngr.start)
    minend = min(rngl.end, rngr.end)
    if not maxstart < minend:
        logger.warning(f"WARNING: no overlap found between {rngl} and {rngr}, double-check input!")

    return impute(cgl.pad(maxstart - rngl.start, rngl.end - minend), cgr.pad(maxstart - rngr.start, rngr.end - minend))



def impute(l: Cigar, r: Cigar):
    """
    This function combines two `Cigar` objects consistisng of the alignment of two queries to a common reference and tries to impute the CIGAR string for a 1:1 alignment of one query to the other.
    By necessity, this is a crude approximation and can in no way replace a full alignment.
    The algorithm assumes the alignments to not start shifted from each other.
    If one of the two alignments has a larger total length than the other one, soft-clipping will be appended to the imputed `Cigar` until the lengths match

    :param l, r: Two `Cigar`s.
    :return: A Cigar containing an approximated alignment of `r` to `l`
    """
    imputed_pairs = []
    if len(l.pairs) == 0 or len(r.pairs) ==0:
        raise ValueError("Empty Cigar input into impute!")

    # prepare the variables for iteration
    liter = iter(l.pairs)
    riter = iter(r.pairs)
    lpr = next(liter)
    rpr = next(riter)
    lprog = 0
    rprog = 0

    while(True):
        try:
            # check if a region has been fully consumed
            if lpr[0]-lprog <= 0 or lpr[1] in cig_clips:
                lprog = 0
                lpr = next(liter)
            if rpr[0]-rprog <= 0 or rpr[1] in cig_clips:
                rprog = 0
                rpr = next(riter)

            # compute the max possible step
            step = min(lpr[0]-lprog, rpr[0]-rprog)
            # compute the type and which strand to step
            t, stepl, stepr = impute_type(lpr[1], rpr[1])
            if t is not None:
                imputed_pairs.append([step, t])
            # do the step #TO/DO make this branchless with cdefs by multiplying with stepl instead
            if stepr:
                rprog += step
            if stepl:
                lprog += step

        except StopIteration:
            break

    # append clippings at end for the sequence that hasn't ran out
    for p in liter:
        imputed_pairs.append([p[0], 'S'])
    for p in riter:
        imputed_pairs.append([p[0], 'S'])

    # return
    cg = Cigar(imputed_pairs)
    cg.clean()
    return cg

cdef impute_type(ltp, rtp):
    """
    Helper function for impute().
    Given two CIGAR types, outputs the one that is imputed
    :param: two valid CIGAR chars (M=IDSHX)
    :return: the CIGAR char imputed for the inputs as well as two bools specifying whether to step ahead on the left/right Cigar.
    """
    # shortcut for the most common part
    if ltp == '=' and rtp == '=':
        return '=', True, True

    if ltp in cig_aln_types:
        if rtp == 'D':
            return 'D', True, True

        if rtp == 'I':
            return 'I', False, True

        if ltp == 'X':
            return 'X' if rtp == '=' else 'M', True, True
        elif ltp == '=': # rtp cannot be '=' in that case b/c that's handled above
            return 'X' if rtp == 'X' else 'M', True, True
        else: # ltp must be 'M'
            return 'M', True, True

    if ltp == 'I':
        if rtp == 'I':
            return 'M', True, True
        else:
            return 'D', True, False # right has deletion relative to left
        
    if ltp == 'D':
        if rtp == 'D':
            return None, True, True # None to skip this region
        elif rtp in cig_aln_types:
            return 'I', True, True
        elif rtp == 'I':
            return 'I', False, True
