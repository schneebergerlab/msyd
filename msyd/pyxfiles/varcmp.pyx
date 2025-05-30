#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

#DEPRECATED

from collections import defaultdict

from msyd.coords import Range
from msyd.multisyn import Multisyn
#from msyd.vars import SNV
import msyd.util as util

logger = util.CustomFormatter.getlogger(__name__)

def get_range_snvs(snvs, rngs):
    """
    A function that finds all SNPs annotated within one of the input region, as a dictionary keyed by organism.
    :params:
        `syri`: dataframe containing the output of `io.extract_syri_snvs`
        `rngs`: a List of `Range` objects to extract snps from.
    :returns:
        A dictionary mapping a `Multisyn` to a dictionary containing the corresponding SNPs identified on it for each organism that `Multisyn` has a position in.
    """
    ret = defaultdict(lambda x: [])
    for _, snv in snvs.iterrows():
        #snv = snv[0]
        if any([snv in rng for rng in rngs]):
            ret[snv.org] += [snv]

    return ret


#TODO implement
# a function getting all snps that are syntenic to a region
def get_syntenic_snps(syri, rng):
    pass

# maybe optionally save with variant annotation later?
class Neighbourhood:
    """
    A class that stores a genomic syntenic neighbourhood, defined by a left and right `Multisyn` neighbour.
    Optionally also stores the syntenic regions overlapping the feature.
    All of `left`, `right` and `overlapping` can be `None`; this means that the neighbourhood is at the leftmost/rightmost edge of a chromosome or that no overlap information is stored, respectively.
    If overlap information is stored, but there are no overlaps, `[]` should be used instead.
    """
    def __init__(self, region:Range, left: Multisyn, right: Multisyn, overlapping):
        self.region = region
        self.left = left
        self.right = right
        self.overlapping = overlapping

    def __eq__(l, r):
        # do not check for overlapping
        return l.left == r.left and l.right == r.right

    ## these two functions compute the gap between the region and the left/right neighbour
    ## 
    def get_left_gap(self):
        # use left/rightmost because region could be inverted
        return self.region.get_leftmost() - self.left.ranges_dict[self.region.org].end

    def get_right_gap(self):
        return self.right.ranges_dict[self.region.org].start - self.region.get_rightmost()

    def get_gaps(self):
        return self.get_left_gap(), self.get_right_gap()

    def dist(l, r):
        # TODO implement a distance function between two neighbourhood objects
        # things to take into account:
        #    – distances to left/right neighbour
        #    – overlapping regions (?)
        # maybe rewrite __eq__ to use dist with a threshold?
        pass

    def get_multisyn_neighbours(rng: Range, multisyns, ref=True, overlapping=True):
        """
        Given a `Range` and a list of multisyntenic regions (as `Multisyn`), returns the regions neighbouring or overlapping with that region on the organism specified in the `Range`.
        :returns: a `Neighbour` object
        """
        # neighbouring := nearest in each direction
        if len(multisyns) < 2:
            logger.warning("get_multisyn_neighbours called with just one multisyn, returning as neighbours!")
            return multisyns

        org = rng.org
        prev = None
        paniter = iter(multisyns.ranges)
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
            logger.warning("get_multisyn_neighbours found rng to be right of all multisyntenic regions! Double-Check if calling with the right set of multisyns")

        return Neighbourhood(rng, left, right, overlapping)
