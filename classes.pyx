#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import copy
import functools
from cigar import Cigar

# these classes form a part of the general SV format
# A position is specified by the organism, chromosome, haplotype and base position
# A range takes a start and an end position. If the end < start, the range is defined as inverted
#TODO use proper cython types, e.g. char for haplo

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
class Position:
    def __init__(self, org:str, chr:int, haplo:str, pos: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.pos = pos

    def __repr__(self):
        return f"Position({self.org}, {self.chr}, {self.haplo}, {self.pos})"

    def __eq__(l, r):
        return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
                l.pos == r.pos

    def __lt__(l, r):
        if l.org != r.org:
            raise ValueError("Comparison between different organisms!")
        if l.chr < r.chr:
            return True
        elif l.chr == r.chr:
            return l.pos < r.pos
        else:
            return False

    def __hash__(self):
        return hash(self.org) + hash(self.chr) + hash(self.haplo) + hash(self.pos)

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
class Range:
    def __init__(self, org:str, chr:int, haplo:str, start: int, end: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.start = start
        self.end = end

    def __repr__(self):
        return f"Range({self.org}, {self.chr}, {self.haplo}, {self.start}, {self.end})"
    
    #def __eq__(l, r):
    #    return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
    #            l.start == r.start & l.start == r.start

    # this operator sorts according to the END, not start value,
    # to enable the end ratchet to work properly
    # TO/DO possible refactor: sort by start here, invert in sorting for algorithm
    # shouldn't really matter as the regions are nonoverlapping, but...
    def __lt__(l, r):
        if l.org != r.org:
            raise ValueError("Comparison between different organisms!")

        if l.chr < r.chr:
            return True
        elif l.chr == r.chr:
            if l.end < r.end:
                return True
            elif l.end == r.end:
                return l.start < r.start
        return False

    def __len__(self):
        return self.end - self.start + 1 # start is inclusive

    def __hash__(self):
        return hash(self.org) + hash(self.chr) + hash(self.haplo) + hash(self.start) + hash(self.end)

    def get_rightmost(self):
        return max(self.start, self.end)

    def get_leftmost(self):
        return min(self.start, self.end)

    def drop(self, start, end):
        """
        :param: 'start'/'end' specify how much to drop on each end.
        :return: A Range missing the specified area.
        """
        if start > len(self) or end > len(self):
            raise ValueError("ERROR: tried to drop more than Range length!")
        return Range(self.org, self.chr, self.haplo, self.start + start, self.end - end)


# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering
class Pansyn:
    """
    A class representing a region syntenic among a set of genomes.
    The parameter `ranges_dict` is a dictionary of genomic `synctools.Range`es storing the location this syntenic region has on each organism.
    This dictionary cannot be None.
    Other attributes are `ref`, which stores the position on the reference -- this attribute can be `None` if using a reference-free algorithm, but none have been implemented so far.
    `cigars_dict` contains a dictionary of `cigar.Cigar` objects corresponding to the alignment of each `Range` to the reference.
    The attribute can be `None` if using approximate position calculation (usually for performance/memory reasons).\n
    Keys in `cigars_dict` correspond with indices in `ranges_dict`.\n
    In the future, `cigars_dict` may also be used for storing pairwise alignments of the core syntenic regions to each other.
    Also, a separate field for an MSA may be added.

    Pansyn implements comparison operators to enable sorting according to the end on the reference.
    For sorting, the `ref` field needs to be set.

    """
    # ranges_dict, cigars_dict have type Dict[String, Range]/Dict[String, Cigar], respectively, but cython cannot deal with generic type hints
    def __init__(self, ref:Range, ranges_dict, cigars_dict):
        #if not ranges_dict:
        #    raise ValueError(f"ERROR: Trying to initialiase Pansyn with no non-reference Range (ref: {ref})")
        #if cigars_dict and not ranges_dict.keys() == cigars_dict.keys():
        #    raise ValueError(f"ERROR: Trying to initialise Pansyn with ranges_dict keys {ranges_dict.keys()} not matching cigars_dict keys {cigars_dict.keys()}!")
        self.ref = ref # optional if using a reference-free algorithm. NONE CURRENTLY IMPLEMENTED!
        self.ranges_dict = ranges_dict
        self.cigars_dict = cigars_dict # optional if using approximate matching

    def __repr__(self):
        return f"Pansyn({self.ref}, {self.ranges_dict})"

    def __eq__(l, r):
        return l.ref == r.ref and l.ranges_dict == r.ranges_dict and l.cigars_dict == r.cigars_dict
        
    # for now, only sorts on the reference (falling back to the Range comparison operator)
    def __lt__(l, r):
        if not l.ref or not r.ref:
            raise ValueError(f"ERROR comparing {l} with {r}: both need to have a reference!")
        return l.ref < r.ref

    def __hash__(self):
        return hash(self.ref)# + hash(self.ranges_dict) + hash(self.cigars_dict) # caused problems with deque

    def add(self, rng:Range, cg: Cigar):
        self.ranges_dict[rng.org] = rng
        if cg:
            if self.cigars_dict:
                self.cigars_dict[rng.org] = cg
            else:
                print("WARNING: attempted to add cigar to Pansyn without cigars_dict, ignoring")

    def get_degree(self):
        return len(self.ranges_dict)

    def get_organisms(self):
        return self.ranges_dict.keys()

    def get_ranges(self):
        return self.ranges_dict.values()

    def get_lens(self):
        return {org: len(self.ranges_dict[org]) for org in self.get_organisms()}

    def check(self):
        """ A function to check a Pansyn object for intactness, mainly for debugging purposes.
        Raises an error if any invariants are violated.
        :returns: None
        """
        if not self.ranges_dict:
            raise ValueError("ERROR in Pansyn.check()! ranges_dict None!")

        if not self.cigars_dict:
            return None

        if self.ranges_dict.keys() != self.cigars_dict.keys():
            raise ValueError("ERROR in Pansyn.check()! ranges_dict keys not matching cigars_dict keys!")
        
        reflen = len(self.ref)
        for org in self.get_organisms():
            if self.cigars_dict[org].get_len(ref=True) != reflen:
                raise ValueError("ERROR in Pansyn.check()! CIGAR length not matching reference length!")
            if self.cigars_dict[org].get_len(ref=False) != len(self.ranges_dict[org]):
                raise ValueError("ERROR in Pansyn.check()! CIGAR length not matching query length!")



    def __add__(self, other):
        """
        Convenience function to concatenate two `Pansyn` objects.
        Uses a shallow copy of the cigar/range to stay without side effects.
        """
        rngs = copy.copy(self.ranges_dict)
        rngs.update(other.ranges_dict)

        cgs = None
        if self.cigars_dict and other.cigars_dict:
            cgs = copy.copy(self.cigars_dict)
            cgs.update(other.cigars_dict)
        elif self.cigars_dict or other.cigars_dict:
            print(f"WARN: Trying to add two Pansyns {self}, {other} with one having CIGARs and one not! Discarding CIGARS!")

        return Pansyn(self.ref, rngs, cgs)

    def drop(self, start, end):
        """
        Returns a new `Pansyn` object with `start`/`end` positions from the start/end of this pansyntenic region removed, respecting cigar alignments if not `None`.
        """
        ref = self.ref.drop(start, end)
        ranges_dict = dict()
        cigars_dict = None
        if not self.cigars_dict:
            ranges_dict = {org:rng.drop(start, end) for (org, rng) in self.ranges_dict.items() if len(rng) > start + end}
        else:
            cigars_dict = dict()
            for org, rng in self.ranges_dict.items():
                cg = self.cigars_dict[org]
                try:
                    start_dropped, cg = cg.get_removed(start, start=True, ref=True)
                    end_dropped, cg = cg.get_removed(end, start=False, ref=True)
                except ValueError:
                    print(f"ERROR: invalid input to cg.get_removed({start}/{end}) on {rng} (len: {len(rng)}). Check if start, end are correct!")
                    continue
                ranges_dict[org] = rng.drop(start_dropped, end_dropped)
                cigars_dict[org] = cg

        return Pansyn(ref, ranges_dict, cigars_dict)

