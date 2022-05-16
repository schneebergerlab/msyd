#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import re

## copied from myUsefulFunctions.py
inttr = lambda x: [int(x[0]), x[1]]
def cgtpl(cg):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    return [inttr(i.split(';')) for i in cg.split(',')[:-1]]
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


reffwd = set(['M', 'D', 'N', '=', 'X'])
qryfwd = set(['M', 'I', 'S', '=', 'X'])
def cg_rem(cg, n, ref=True, start=True):
    """
    Takes cigar as input, removes from 'ref'erence/query strand until 'n' bases from the OTHER strand have been removed.
    Starts at the 'start'/end.
    :return: the cigar string with these bases removed
    """
    #if type(cg) == str:
    #    cg = cgtpl(cg)

    ind = 0 if start else -1
    skip = 0
    fwd = reffwd if ref else qryfwd
    altfwd = qryfwd if ref else reffwd

    while n > 0:
        cgi = cg[ind]
        print(n, start, cgi)
        if cgi[1] in altfwd:
            skip += cgi[0]

        if cgi[1] not in fwd:
            if start:
                cg = cg[1:]
            else:
                cg = cg[:-1]
            continue

        # cgi must be in fwd
        n -= cgi[0]
        if n >= 0:
            if start:
                cg = cg[1:]
            else:
                cg = cg[:-1]
        else:
            cgi[0] = -n # n == n- cg[ind][0] implies -n == cg[ind][0] -n
            if cgi[1] in altfwd: # subtract the overcounting
                skip += n

    return (skip, cg)


if __name__ == "__main__":
    cigar = "10=5D3X7=4I8="
    import sys
    print(cg_rem(cigar, int(sys.argv[1])))
    print(cg_rem(cigar, int(sys.argv[1]), ref=False))
    print(cg_rem(cigar, int(sys.argv[1]), start=False))
