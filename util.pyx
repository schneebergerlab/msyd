#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import re
import copy
import multiprocessing
import functools

## constants
reffwd = set(['M', 'D', 'N', '=', 'X'])
qryfwd = set(['M', 'I', 'S', '=', 'X'])
inttr = lambda x: [int(x[0]), x[1]]

class Cigar:
    """
    A CIGAR string represented as a list of lists.
    """

    def __init__(self, pairs):
        self.pairs = pairs

    ## copied from myUsefulFunctions.py
    def from_string(cg):
        """
        Takes a cigar string as input and returns a cigar tuple
        """
        
        for i in "MIDNSHPX=":
            cg = cg.replace(i, ';'+i+',')
        return Cigar([inttr(i.split(';')) for i in cg.split(',')[:-1]])

    def cg_genlen(self, ref=True):
        """
        Takes cigar as input, and return the number of bases covered by it the reference or query genome.
        """
        s = reffwd if ref else qryfwd
        return sum([int(i[0]) for i in cg if i[1] in s])

    def get_removed(self, n, ref=True, start=True):
        """
        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
        """
        cg = copy.deepcopy(self) #TODO possibly very slow

        ind = 0 if start else -1
        skip = 0
        fwd = reffwd if ref else qryfwd
        altfwd = qryfwd if ref else reffwd

        while n > 0:
            cgi = cg.pairs[ind]
            #print(n, start, cgi)
            if cgi[1] in altfwd:
                skip += cgi[0]

            if cgi[1] not in fwd:
                # if ever removing the copy, refactor!
                del cg.pairs[ind]
                continue

            # cgi must be in fwd
            n -= cgi[0]
            if n >= 0:
                del cg.pairs[ind]
            else:
                cgi[0] = -n # n == n- cg[ind][0] implies -n == cg[ind][0] -n
                if cgi[1] in altfwd: # subtract the overcounting
                    skip += n

        return (skip, cg)

    def __repr__(self):
        return f"Cigar({self.pairs})"

# copied from https://stackoverflow.com/questions/50878960/parallelize-pythons-reduce-command
def parallel_reduce(reduceFunc, l, numCPUs):
    if numCPUs == 1 or len(l) <= 100:
            returnVal= functools.reduce(reduceFunc, l[1:], l[0])
            return returnVal

    parent1, child1 = multiprocessing.Pipe()
    parent2, child2 = multiprocessing.Pipe()
    p1 = multiprocessing.Process(target=parallel_reduce, args=(reduceFunc, l[:len(l) // 2], numCPUs // 2, child1, ) )
    p2 = multiprocessing.Process(target=parallel_reduce, args=(reduceFunc, l[len(l) // 2:], numCPUs // 2 + numCPUs%2, child2, ) )
    p1.start()
    p2.start()
    leftReturn, rightReturn = parent1.recv(), parent2.recv()
    p1.join()
    p2.join()
    returnVal = reduceFunc(leftReturn, rightReturn)
    return returnVal


if __name__ == "__main__":
    cg = Cigar.from_string("10=5D3X7=4I8=")
    import sys
    print(cg.get_removed(int(sys.argv[1])))
    print(cg.get_removed(int(sys.argv[1]), ref=False))
    print(cg.get_removed(int(sys.argv[1]), start=False))
