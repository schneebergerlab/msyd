#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import copy
import functools
#import cython
import logging

from msyd.cigar import Cigar
import msyd.scripts.util as util

logger = util.CustomFormatter.getlogger(__name__)

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
        #return f"Position({self.org}, {self.chr}, {self.haplo}, {self.pos})"
        return f"Position({self.org}, {self.chr}, {self.pos})"

    def to_pff(self):
        """Transform this `Position` into pansynteny file format
        """
        #return f"{self.chr}:{self.haplo}:{self.pos}"
        return f"{self.chr}:{self.pos}"

    def __eq__(l, r):
        if not isinstance(r, Position):
            return False
        return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
                l.pos == r.pos

    def __lt__(l, r):
        if l.org != r.org:
            logger.error(f"Comparison between different organisms: {l.org} != {r.org}")
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
#@cython.total_ordering
class Range:
    def __init__(self, org:str, chr:str, haplo:str, start: int, end: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.start = start # inclusive
        self.end = end # inclusive

    def __repr__(self):
        #return f"Range({self.org}, {self.chr}, {self.haplo}, {self.start}, {self.end})"
        return f"Range({self.org}, {self.chr}, {self.start}, {self.end})"

    def to_pff(self):
        """Transform this `Range` into the form specified by PFF
        """
        #return f"{self.chr}:{self.haplo}:{self.start}-{self.end}"
        return f"{self.chr}:{self.start}-{self.end}"

    def to_pff_org(self):
        """Transform this `Range` into the form specified by PFF, with the sample name being prepended as specified for the realigned reference haplotype
        """
        #return f"{self.org}:{self.chr}:{self.haplo}:{self.start}-{self.end}"
        return f"{self.org}:{self.chr}:{self.start}-{self.end}"


    def read_pff(org:str, cell: str):
        """Parse a Range in PFF format
        PFF format is :-separated and must contain the chromosome first, then followed by a haplotype (optional) and a start and end position separated by -.
        If the end position is before the start position, the range is treated as inverted.
        `org` specifies the organism the returned range should have. If this range is just used for filtering, None may be passed.
        Examples: Chr1:mat:1000-2000, Chr3:10000-50000
        """
        #TODO error handling in here
        #print(cell)
        cellarr = cell.split(':')
        if len(cellarr) < 2 or len(cellarr) > 4:
            raise ValueError(f"Invalid PFF Range string: {cell}")
        if len(cellarr) == 4: # if a ref name is specified in the cell, that overrides the argument
            org = cellarr[0]
            cellarr = cellarr[1:]
        start = int(cellarr[-1].split('-')[0])
        end = int(cellarr[-1].split('-')[1])
        hapl = cellarr[1] if len(cellarr) == 3 else None
        chrom = cellarr[0] #util.chrom_to_int(cellarr[0])
        # should chr be int'ed as well?
        return Range(org, chrom, hapl, start, end)

    
    def __eq__(l, r):
        if not isinstance(r, Range):
            return False
        return l.org == r.org and l.chr == r.chr and l.haplo == r.haplo and \
                l.start == r.start & l.start == r.start

    # this operator sorts according to the END, not start value,
    # to enable the end ratchet to work properly
    # TO/DO possible refactor: sort by start here, invert in sorting for algorithm
    # shouldn't really matter as the regions are nonoverlapping, but...
    def __lt__(l, r):
        if l.org != r.org:
            logger.error(f"Comparison between different organisms: {l.org} != {r.org}")
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

    def __contains__(self, item):
        #TODO expand to work with variation superclass, how to handle transpositions?
        """Given a `Position` or `Range` object, outputs `True` if the position/region is fully contained by this `Range`. Ignores haplotype information.
            :param item: a genomic feature that is checked for being inside this `Range`.
            :type item: `Range` or `Position`
            :returns: `True` if the supplied feature is within this `Range`.
            :rtype: `Bool`
        """
        if isinstance(item, Position):
            return self.chr == item.chr and self.start <= item.pos <= self.end
        elif isinstance(item, Range):
            return self.chr == item.chr and self.start <= item.start and item.end <= self.end
        else:
            return False

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
        if start < 0 or end < 0:
            raise ValueError("ERROR: tried to drop negative value!")
        return Range(self.org, self.chr, self.haplo, self.start + start, self.end - end)

    def is_inverted(self):
        return self.start <= self.end

    def check(self):
        if self.start < 0:
            return False
        if self.end < 0:
            return False


# TODO: Maybe, it is possible to have the pansyn class in pansyn.pyx or a separate classes.pyx file? Alternatively, it seems that the coords.pyx just defines Class:pansyn. Cant we just move everything to pansyn.pyx
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
        return f"Pansyn({self.ref}, {self.ranges_dict})"#, {self.cigars_dict})"

    def __eq__(l, r):
        if not isinstance(r, Pansyn):
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
                logger.warning("attempted to add cigar to Pansyn without cigars_dict, ignoring")

    def get_degree(self):
        return len(self.ranges_dict)

    def get_orgs(self):
        return self.ranges_dict.keys()
    def get_organisms(self):
        return self.ranges_dict.keys()

    def get_ranges(self):
        return self.ranges_dict.values()

    def get_lens(self):
        return {org: len(self.ranges_dict[org]) for org in self.get_organisms()}

    def check(self):
        """
        A function to check a Pansyn object for intactness, mainly for debugging purposes.
        Returns `False` if any invariants are violated.
        :returns: `True` if the object is a valid `Pansyn` object, else `False`
        """

        if not self.ranges_dict:
            return False
            #raise ValueError("ERROR in Pansyn.check()! ranges_dict None!")

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
            #raise ValueError("ERROR in Pansyn.check()! ranges_dict keys not matching cigars_dict keys!")
        
        # length check
        reflen = len(self.ref)
        for org in self.get_organisms():
            if self.cigars_dict[org].get_len(ref=True) != reflen:
                return False
                #raise ValueError("ERROR in Pansyn.check()! CIGAR length not matching reference length!")
            if self.cigars_dict[org].get_len(ref=False) != len(self.ranges_dict[org]):
                return False
                #raise ValueError("ERROR in Pansyn.check()! CIGAR length not matching query length!")


    def __add__(self, other):
        """
        Convenience function to concatenate two `Pansyn` objects.
        Uses a shallow copy of the cigar/range to stay without side effects.
        """
        # if this is or other is an empty pansyn, return early
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
            logger.warning(f"Trying to add two Pansyns {self}, {other} with one having CIGARs and one not! Discarding CIGARS!")

        return Pansyn(self.ref, rngs, cgs)

    def drop(self, start, end):
        #, prop=False):
        #:param prop: Controls whether to drop the same amount in absolute terms (default) or proportional to the region lengths when dropping from a `Pansyn` without CIGAR strings.
        """
        Returns a new `Pansyn` object with `start`/`end` positions from the start/end of this pansyntenic region removed, respecting cigar alignments if not `None`.
        """
        if start < 0 or end < 0 or start + end > len(self.ref):
            logger.error(f"Tried to drop invalid start ({start}) or end ({end}) on this Pansyn with length on the reference {len(self.ref)}")
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

                if start + end < len(rng): # TODO maybe handle this case to drop proportionally, i.e. if the drop is 10% of ref, drop 10% of qry instead of having absolute length the same
                    ranges_dict[org] = rng.drop(start, end)
        else:
            cigars_dict = dict()
            for org, rng in self.ranges_dict.items():
                cg = self.cigars_dict[org]
                try:
                    if not cg.is_empty():
                        start_dropped, cg = cg.get_removed(start, start=True, ref=True)
                    else:
                        logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}). Skipping!")
                        continue
                    if not cg.is_empty():
                        end_dropped, cg = cg.get_removed(end, start=False, ref=True)
                    else:
                        logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}). Skipping!")
                        continue
                except ValueError:
                    logger.warning(f"Tried to drop more({start}/{end}) than length on {rng}(len: {len(rng)}) on org {org}. Skipping!")
                    continue
                ranges_dict[org] = rng.drop(start_dropped, end_dropped)
                cigars_dict[org] = cg

        return Pansyn(ref, ranges_dict, cigars_dict)

