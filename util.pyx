#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

## copied from myUsefulFunctions.py
def cgtpl(cg):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    return [i.split(';') for i in cg.split(',')[:-1]]
#end

def cggenlen(cg, gen):
    """
    Takes cigar as input, and return the number of bases covered by it the reference
    or query genome.
    Cigar strings are first converted to cigar tuples.
    """
    if type(cg) == str:
        cg = cgtpl(cg)
    if gen not in ['r', 'q']:
        raise ValueError('gen need to "r" or "q" for reference or query')
        return
    s = set(['M', 'D', 'N', '=', 'X']) if gen == 'r' else set(['M', 'I', 'S', '=', 'X'])
    l = sum([int(i[0]) for i in cg if i[1] in s])
    return l
#end


refforw = set(['M', 'D', 'N', '=', 'X'])
qryforw = set(['M', 'I', 'S', '=', 'X'])
def remfromcg(cg, n, ref=True, start=True):
    """
    Takes cigar as input, removes from 'ref'erence/query strand until 'n' bases from the OTHER strand have been removed.
    Starts at the 'start'/end.
    :return: the cigar string with these bases removed
    """
    if type(cg) == str:
        cg = cgtpl(cg)

    if not start:
        cg = cg[::-1]

    while n > 0:
        pass

