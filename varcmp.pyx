#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import numpy as np
from pansyn import Range, Pansyn
import pandas as pd

def get_snps(syri, rngs):
    """
    A function that finds all SNPs annotated within one of the input region, as a dictionary keyed by organism.
    :params:
        `syri`: dataframe containing the information from a syri file, as parsed by `ingest.readsyriout`
        `rngs`: a List of `Range` objects to extract snps from.
    :returns:
        A dictionary mapping a `Pansyn` to a dictionary containing the corresponding SNPs identified on it for each organism that `Pansyn` has a position in.
    """
    pass

# maybe optionally save with variant annotation later?
class Neighbourhood:
    """
    A class that stores a genomic syntenic neighbourhood, defined by a left and right `Pansyn` neighbour.
    Optionally also stores the syntenic regions overlapping the feature.
    All of `left`, `right` and `overlapping` can be `None`; this means that the neighbourhood is at the leftmost/rightmost edge of a chromosome or that no overlap information is stored, respectively.
    If overlap information is stored, but there are no overlaps, `[]` should be used instead.
    """
    def __init__(self, left: Pansyn, right: Pansyn, overlapping):
        self.left = left
        self.right = right
        self.overlapping = overlapping

    def __eq__(l, r):
        # do not check for overlapping
        return l.left == r.left and l.right == r.right

    def dist(l, r):
        # TODO implement a distance function between two neighbourhood objects
        # things to take into account:
        #    – distances to left/right neighbour
        #    – overlapping regions (?)
        # maybe rewrite __eq__ to use dist with a threshold?
        pass

    def get_pansyn_neighbours(rng: Range, pansyns, ref=True, overlapping=True):
        """
        Given a `Range` and a list of pansyntenic regions (as `Pansyn`), returns the regions neighbouring or overlapping with that region on the organism specified in the `Range`.
        :returns: a `Neighbour` object
        """
        # neighbouring := nearest in each direction
        if len(pansyns) < 2:
            print("WARN: get_pansyn_neighbours called with just one pansyn, returning as neighbours!")
            return pansyns

        org = rng.org
        prev = None
        paniter = iter(pansyns.ranges)
        prev = next(paniter)
        cur = next(paniter)
        left = None
        right = None
        overlapping = [] if overlapping else None
        try:
            while cur.ref.end < rng.start if ref else cur.ranges_dict[org].end < rng.start: # maybe TODO refactor to use more elegant polymorphism than ternaries?
                prev = cur
                cur = next(paniter)

            left = prev
            del prev # no longer need to save previous values
            
            # get overlapping 'neighbours'
            while cur.ref.start < rng.end if ref else cur.ranges_dict[org].start < rng.end:
                if overlapping:
                    overlapping.append(cur)
                cur = next(paniter)

            right = cur
                    
        except StopIteration:
            print("WARN: get_pansyn_neighbours found rng to be right of all pansyntenic regions! Double-Check if calling with the right set of pansyns")

        return Neighbourhood(left, right, overlapping)
