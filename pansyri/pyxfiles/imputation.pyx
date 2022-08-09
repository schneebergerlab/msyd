#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

from pansyn.classes.cigar import Cigar, cig_clips, cig_aln_types

"""
TODOs
– write unit tests
– write frontend:
    • impute including ranges, to handle different positions (use clipping)
    • impute unequal sequences, also use clipping
    • maybe have a convenience function imputing everything given two bam files?
    • maybe have a convenience function automatically detecting core/cross synteny from such imputed bams?
    • maybe have a fucntion that automatically imputes along all pansyntenic regions?
    => convenience functions may be moved to util
"""


def impute(l, r):
    """
    This function combines two CIGAR strings of two queries on one reference and tries to impute the CIGAR string for a 1:1 alignment of one query to the other.
    By necessity, this is a crude approximation and can in no way replace a full alignment.
    The algorithm assumes the alignments to not start shifted from each other.
    :param: Two Cigars.
    :return: A Cigar.
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
            t, stepl, stepr = Cigar.impute_type(lpr[1], rpr[1])
            if t is not None:
                imputed_pairs.append([step, t])
            # do the step #TO/DO make this branchless with cdefs
            if stepr:
                rprog += step
            if stepl:
                lprog += step

        except StopIteration:
            break

    # append clippings at end for the sequence that hasn't ran out
    #for p in liter:
    #    imputed_pairs.append([p[0], 'S'])
    #for p in riter:
    #    imputed_pairs.append([p[0], 'S'])

    # return
    cg = Cigar(imputed_pairs)
    cg.clean()
    return cg

def impute_type(ltp, rtp):
    """
    Helper function for impute().
    Given two CIGAR types, outputs the one that is imputed
    :param: two valid CIGAR chars (M=IDSHX)
    :return: the CIGAR char imputed for the inputs as well as two bools specifying whether to step ahead on the left/right Cigar.
    """
    if ltp in cig_aln_types:
        if rtp == 'D':
            return 'D', True, True

        if rtp == 'I':
            return 'I', True, True

        if ltp == 'X' or rtp == 'X':
            return 'X', True, True
        else:
            if ltp == '=' and rtp == '=':
                return '=', True, True
            else: # at least one must be an M, none is an X
                return 'M', True, True

    if ltp == 'I':
        if rtp == 'I': # TO/DO be even more conservative and resolve this to D followed by I?
            return 'M', True, True
        else:
            return 'D', True, False # right has deletion relative to left
        
    if ltp == 'D':
        if rtp == 'D':
            return None, True, True # None to skip this region
        elif rtp in cig_aln_types:
            return 'I', True, True
        elif rtp == 'I':
            return 'I', True, True
