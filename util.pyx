#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
import multiprocessing

# copied from https://stackoverflow.com/questions/50878960/parallelize-pythons-reduce-command
# doesn't seem to be very fast?
def parallel_reduce(reduceFunc, l, numCPUs):
    if numCPUs == 1 or len(l) <= 100:
            returnVal = functools.reduce(reduceFunc, l[1:], l[0])
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


