#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import copy
import functools
import traceback
#import cython
import logging

from libcpp.vector cimport vector

from msyd.cigar import Cigar
import msyd.util as util
from msyd.coords import Range

logger = util.CustomFormatter.getlogger(__name__)

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering
#cdef
class Multisyn:
    #    cdef:
    #        Range ref
    #        dict ranges_dict
    #        dict cigars_dict
    """
    A class representing a region syntenic among a set of genomes.
    The parameter `ranges_dict` is a dictionary of genomic `synctools.Range`es storing the location this syntenic region has on each organism.
    This dictionary cannot be None.
    Other attributes are `ref`, which stores the position on the genome used as a reference in calling this Multisyn. Different multisyn may have different references, and coordinates are only comparable between the same reference (as stored in `ref.org`).
    `cigars_dict` contains a dictionary of `cigar.Cigar` objects corresponding to the alignment of each `Range` to the reference.
    The attribute can be `None` if using approximate position calculation (usually for performance/memory reasons).\n
    Keys in `cigars_dict` correspond with indices in `ranges_dict`.\n
    In the future, `cigars_dict` may also be used for storing pairwise alignments of the core syntenic regions to each other.
    Also, a separate field for an MSA may be added.

    Multisyn implements comparison operators to enable sorting according to the end on the reference.
    For sorting, the `ref` field needs to be set.

    """


    # ranges_dict, cigars_dict have type Dict[String, Range]/Dict[String, Cigar], respectively, but cython cannot deal with generic type hints
    def __init__(self, ref:Range, ranges_dict, cigars_dict):
        #if not ranges_dict:
        #    raise ValueError(f"ERROR: Trying to initialiase Multisyn with no non-reference Range (ref: {ref})")
        #if cigars_dict and not ranges_dict.keys() == cigars_dict.keys():
        #    raise ValueError(f"ERROR: Trying to initialise Multisyn with ranges_dict keys {ranges_dict.keys()} not matching cigars_dict keys {cigars_dict.keys()}!")
        self.ref = ref # optional if using a reference-free algorithm. NONE CURRENTLY IMPLEMENTED!
        self.ranges_dict = ranges_dict
        self.cigars_dict = cigars_dict # optional if using approximate matching

    def __repr__(self):
        return f"Multisyn({self.ref}, {self.ranges_dict})"#, {self.cigars_dict})"

    def __eq__(l, r):
        if not isinstance(r, Multisyn):
            return False
        return l.ref == r.ref and l.ranges_dict == r.ranges_dict and l.cigars_dict == r.cigars_dict
        
    # for now, only sorts on the reference (falling back to the Range comparison operator)
    def __lt__(l, r):
        if not l.ref or not r.ref:
            logger.error(f"comparing {l} with {r}: both need to have a reference!")
            raise ValueError(f"ERROR comparing {l} with {r}: both need to have a reference!")
        if l.ref.org != r.ref.org:
            logger.error(f"Comparison between different references trying to compare {l} and {r}")
            raise ValueError("Comparison between different references!")
        return l.ref < r.ref

    def __hash__(self):
        return hash(self.ref)# + hash(self.ranges_dict) + hash(self.cigars_dict) # caused problems with deque

    def add(self, rng:Range, cg: Cigar):
        self.ranges_dict[rng.org] = rng
        if cg:
            if self.cigars_dict:
                self.cigars_dict[rng.org] = cg
            else:
                logger.warning("attempted to add cigar to Multisyn without cigars_dict, ignoring")

    def get_degree(self):
        return len(self.ranges_dict) + 1 # count the ref as well

    def get_orgs(self):
        return self.get_organisms()
    def get_organisms(self):
        return self.ranges_dict.keys()

    def get_ranges(self):
        return self.ranges_dict.values()

    def get_lens(self):
        return {org: len(self.ranges_dict[org]) for org in self.get_organisms()}

    def check(self):
        """
        A function to check a Multisyn object for intactness, mainly for debugging purposes.
        Returns `False` if any invariants are violated.
        :returns: `True` if the object is a valid `Multisyn` object, else `False`
        """

        if not self.ranges_dict:
            return False
            #raise ValueError("ERROR in Multisyn.check()! ranges_dict None!")

        if not self.ref.check() or self.ref.is_inverted():
            return False

        for rng in self.ranges_dict.values():
            if not rng.check() or rng.is_inverted():
                return False


        ## CIGAR checks
        if not self.cigars_dict:
            return True

        if self.ranges_dict.keys() != self.cigars_dict.keys():
            return False
            #raise ValueError("ERROR in Multisyn.check()! ranges_dict keys not matching cigars_dict keys!")
        
        # length check
        reflen = len(self.ref)
        for org in self.get_organisms():
            if self.cigars_dict[org].get_len(ref=True) != reflen:
                return False
                #raise ValueError("ERROR in Multisyn.check()! CIGAR length not matching reference length!")
            if self.cigars_dict[org].get_len(ref=False) != len(self.ranges_dict[org]):
                return False
                #raise ValueError("ERROR in Multisyn.check()! CIGAR length not matching query length!")

    def trim_matching_inplace(self):
        """
        Trims the alignment of this multisyn to remove mismatches at the start/end.
        Mutates self!
        """
        if not self.cigars_dict: # can't trim if there are no alignments
            return

        if len(self.ref) > 2000:
            print(self)
            print(self.cigars_dict)

        trims = [cg.trim_matching(only_pos=True) for cg in self.cigars_dict.values()]
        start = max(tup[2] for tup in trims)
        end = max(tup[3] for tup in trims)
        if start > 0 or end > 0:
            if len(self.ref) > 2000:
                print(start, end)
            self.drop_inplace(start, end)


    def split_indels(self, thresh):
        """

        """
        if len(self.ranges_dict) > 1:
            logger.error("Splitting multisyns with more than 2 genomes not supported! Split before intersecting")
            return ValueError("Splitting multisyns with more than 2 genomes not supported! Split before intersecting")

        org = list(self.ranges_dict)[0]
        altrng = self.ranges_dict[org]
        splits = list(self.cigars_dict.values())[0].split_indels(thresh=thresh)

        # shortcut to avoid copying unnecessarily
        if len(splits) == 0:
            return [self]

        return [Multisyn(
                        Range(self.ref.org, self.ref.chr, split[0], split[1]),
                        {org: Range(altrng.org, altrng.chr, split[2], split[3])},
                        {org: split[4]}
                        )
                    for split in splits]

    def __add__(self, other):
        """
        Convenience function to concatenate two `Multisyn` objects.
        Uses a shallow copy of the cigar/range to stay without side effects.
        """
        # if this is or other is an empty multisyn, return early
        if not other.ranges_dict:
            return self
        elif not self.ranges_dict:
            ret = copy.copy(other)
            ret.ref = self.ref # not sure if necessary, but doesn't hurt I guess
            return ret

        rngs = copy.copy(self.ranges_dict)
        rngs.update(other.ranges_dict)

        cgs = None
        if self.cigars_dict and other.cigars_dict:
            cgs = copy.copy(self.cigars_dict)
            cgs.update(other.cigars_dict)
        elif self.cigars_dict or other.cigars_dict:
            logger.warning(f"Trying to add two Multisyns {self}, {other} with one having CIGARs and one not! Discarding CIGARS!")

        return Multisyn(self.ref, rngs, cgs)


    def drop_on_org(self, start, end, org):
        """
        Similar to drop, but supports specifying the number of bases to drop on an organism other than the reference.
        Returns a new `Multisyn` object with `start`/`end` positions from the start/end of this multisyntenic region removed counting on `org`, respecting cigar alignments if not `None`.
        If this `Multisyn` has no cigar object, it will be identical to a drop on the reference (as the genomes are assumed to align perfectly).
        """
        # in case its called with the ref org already, directly drop
        # this fn shouldn't usually be called like this, but handle anyway
        if org == self.ref.org:
            return self.drop(start, end)

        assert(org in self.ranges_dict)
        # get no of bases corresponding to this drop on the reference
        if self.cigars_dict:
            start, end = self.cigars_dict[org].trim(start, end, only_pos=True, ref=False)

        return self.drop(start, end)

    def drop_on_org_inplace(self, start, end, org):
        """
        Same as `drop_on_org`, but instead of returning a new Multisyn object, this mutates the current object.
        """
        # in case its called with the ref org already, directly drop
        # this fn shouldn't usually be called like this, but handle anyway
        if org == self.ref.org:
            self.drop_inplace(start, end)
            return

        assert(org in self.ranges_dict)
        # get no of bases corresponding to this drop on the reference
        if self.cigars_dict:
            start, end = self.cigars_dict[org].trim(start, end, only_pos=True, ref=False)

        self.drop_inplace(start, end)

    def drop(self, start, end):
        #, prop=False): #DEPRECATED
        #:param prop: Controls whether to drop the same amount in absolute terms (default) or proportional to the region lengths when dropping from a `Multisyn` without CIGAR strings.
        """
        Returns a new `Multisyn` object with `start`/`end` positions from the start/end of this multisyntenic region removed, respecting cigar alignments if not `None`.
        """
        if start < 0 or end < 0 or start + end > len(self.ref):
            logger.error(f"Tried to drop invalid start ({start}) or end ({end}) on this Multisyn with length on the reference {len(self.ref)}")
            raise ValueError("tried to drop invalid start/end!")

        ref = self.ref.drop(start, end)
        #reflen = len(ref)
        #pstart = start/reflen
        #pend = end/reflen

        ranges_dict = dict()
        cigars_dict = None
        if not self.cigars_dict:
            for org, rng in self.ranges_dict.items():
                #if prop:
                #    l = len(rng)
                #    start = int(pstart*l)
                #    end = int(pend*l)

                if start + end < len(rng):
                    ranges_dict[org] = rng.drop(start, end)
        else:
            cigars_dict = dict()
            for org, rng in self.ranges_dict.items():
                cg = self.cigars_dict[org]
                try:
                    if not cg.is_empty():
                        start_dropped, cg = cg.get_removed(start, start=True, ref=True)
                    else:
                        print(traceback.format_exc())
                        logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}). Skipping!")
                        continue
                    if not cg.is_empty():
                        end_dropped, cg = cg.get_removed(end, start=False, ref=True)
                    else:
                        print(traceback.format_exc())
                        logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}). Skipping!")
                        continue

                except ValueError:
                    print(traceback.format_exc())
                    logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}) on org {org}. Skipping!")
                    continue

                ranges_dict[org] = rng.drop(start_dropped, end_dropped)
                cigars_dict[org] = cg

        return Multisyn(ref, ranges_dict, cigars_dict)

    #TODO? write test testing that this is equivalent to drop
    def drop_inplace(self, start, end):
        #, prop=False): #DEPRECATED
        #:param prop: Controls whether to drop the same amount in absolute terms (default) or proportional to the region lengths when dropping from a `Multisyn` without CIGAR strings.
        """
        Performs the same function as `drop`, but mutates this object instead of returning a new one.
        Mutates this `Multisyn` object to remove `start`/`end` positions from the start/end, respecting cigar alignments if not `None`.
        """
        if start < 0 or end < 0 or start + end > len(self.ref):
            logger.error(f"Tried to drop invalid start ({start}) or end ({end}) on this Multisyn with length on the reference {len(self.ref)}")
            raise ValueError("tried to drop invalid start/end!")

        self.ref = self.ref.drop(start, end)
        #reflen = len(ref)
        #pstart = start/reflen
        #pend = end/reflen

        if not self.cigars_dict:
            for org, rng in self.ranges_dict.items():
                #if prop:
                #    l = len(rng)
                #    start = int(pstart*l)
                #    end = int(pend*l)

                if start + end < len(rng):
                    self.ranges_dict[org] = rng.drop(start, end)
        else:
            for org, rng in self.ranges_dict.items():
                cg = self.cigars_dict[org]
                try:
                    if not cg.is_empty():
                        start_dropped, cg = cg.get_removed(start, start=True, ref=True)
                    else:
                        print(traceback.format_exc())
                        logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}). Skipping!")
                        continue
                    if not cg.is_empty():
                        end_dropped, cg = cg.get_removed(end, start=False, ref=True)
                    else:
                        print(traceback.format_exc())
                        logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}). Skipping!")
                        continue

                except ValueError:
                    print(traceback.format_exc())
                    logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}) on org {org}. Skipping!")
                    continue

                self.ranges_dict[org] = rng.drop(start_dropped, end_dropped)
                self.cigars_dict[org] = cg



cdef class Multisyn_container:
    """
    A datastructure for storing Multisyn objects on the same chromosome.
    Multisyns are segregated by chromosome to allow efficient lock-free parallelism between different chromosomes.
    Internally, Multisyns are stored in a C++ vector sorted by position on the first chosen reference. Merasyns that do not have a position on the first reference are sorted immediately after the coresyn before them, sorted alphabetically by their reference name and by position on that reference.
    This allows the use of binary search to efficiently retrieve even sequences not present in the first reference.
    This sorting is also used for iterating over the Multisyns during the synteny intersection step.
    """
    cdef list[Multisyn] multisyns

    def __cinit__(self, cap: int):
        """
        `cap` specifies the initial capacity the vector backing should be initialised with.
        """
        self.multisyns = list()
        if cap:
            self.reserve(cap)

    def __len__(self):
        return self.multisyns.size()

    def __repr__(self):
        return self.to_string(10)

    def to_string(self, n: int):
        pass

    def __getitem__(self, key: int):
        """
        `key` may be an index, or a tuple for a slice.
        """
        #TODO handle slicing – copy to smaller Multisyn_container maybe?
        if len(key) > 1:
            raise NotImplemented("Slicing not yet supported!")
        return self.getat(key)

    cdef getat(self, pos:int):
        return self.vector.at(pos)

    cdef copy_slice(self, start:int, end:int):
        pass

    cpdef push_back(self, ms: Multisyn):
        """
        Append a new Multisyn.
        """
        #TODO check if sorted?
        self.vector.push_back(ms)

    cpdef reserve(self, cap:int):
        """
        Pre-allocate to contain `cap` elements.
        """
        self.vector.reserve(cap)

    cpdef shrink_to_fit(self):
        """
        Shrink the backing vector to the current contents.
        Frees up memory, especially if this Container has been overallocated previously.
        """
        self.vector.shrink_to_fit()


    def find(self, org: str, start: int, end: int):
        """
        Return a slice of all multisyns overlapping with `start:end` (inclusive) on `org`.
        """
        #TODO check chromosome?
        endind = self.find_ind(org, end)
        # find the last position still covered here
        while self[endind].ranges_dict[org].end < endind and endind < len(self):
            endind += 1

        #TODO does this work if start is not covered?
        return self[self.find_ind(org, start):endind]

    def find(self, org:str, pos: int):
        """
        Find the Multisyn covering `pos` on `org`.
        Return None if that position is not covered.
        """
        #TODO check chromosome?
        ms = self[self.find_ind(org, pos)]
        return ms if pos in ms else None

    cdef find_ind(self, org:str, pos: int):
        """
        Find the index at or immediately before `pos` using binary search.
        """
        pass


