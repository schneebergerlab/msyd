#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
#import cython
import logging

import msyd.util as util

logger = util.CustomFormatter.getlogger(__name__)

# these classes form a part of the general SV format
# A position is specified by the organism, chromosome, haplotype and base position
# A range takes a start and an end position. If the end < start, the range is defined as inverted

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
cdef class Position:
    cdef:
        public str org
        public str chr
        public unsigned long pos

    def __cinit__(self, org:str = None, chr:str = None, pos:int = 0):
    #def __init__(self, org:str, chr:str, pos:int):
        """
        All args are optional to support pickling; always set them otherwise!
        """
        self.org = org
        self.chr = chr
        self.pos = pos

    def __repr__(self):
        #return f"Position({self.org}, {self.chr}, {self.pos})"
        return f"Position({self.org}, {self.chr}, {self.pos})"

    def to_psf(self):
        """Transform this `Position` into population synteny file format
        """
        return f"{self.chr}:{self.pos}"

    def to_psf_org(self):
        """Transform this `Position` into population synteny file format, including the org
        """
        return f"{self.org}:{self.chr}:{self.pos}"

    # support pickling, for use with multiprocessing
    def __getstate__(self):
        return self.to_psf_org()

    def __setstate__(self, state):
        # kind of ugly, but pickling requires to keep the object
        rng = read_psf_pos(None, state)
        self.org = rng.org
        self.chr = rng.chr
        self.start = rng.pos

    def __eq__(l, r):
        if not isinstance(r, Position):
            return False
        return l.org == r.org and l.chr == r.chr and \
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
        return hash(self.org) + hash(self.chr) + hash(self.pos)

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
#@cython.total_ordering
cdef class Range:
    cdef:
        public str org
        public str chr
        public unsigned long start
        public unsigned long end

    def __cinit__(self, org:str = None, chr:str = None, start:int = 0, end:int = 0):
    #def __init__(self, org:str, chr:str, start:int, end:int):
        """
        All args are optional to support pickling; always set them otherwise!
        """
        self.org = org
        self.chr = chr
        self.start = start # inclusive
        self.end = end # inclusive

    def __repr__(self):
        #return f"Range({self.org}, {self.chr}, {self.start}, {self.end})"
        return f"Range({self.org}, {self.chr}, {self.start}, {self.end})"

    def to_psf(self):
        """Transform this `Range` into the form specified by PSF
        """
        #return f"{self.chr}:{self.start}-{self.end}"
        return f"{self.chr}:{self.start}-{self.end}"

    def to_psf_org(self):
        """Transform this `Range` into the form specified by PSF, with the sample name being prepended as specified for the realigned reference haplotype
        """
        #return f"{self.org}:{self.chr}:{self.start}-{self.end}"
        return f"{self.org}:{self.chr}:{self.start}-{self.end}"

    # support pickling, for use with multiprocessing
    def __getstate__(self):
        return self.to_psf_org()

    def __setstate__(self, state):
        # kind of ugly, but pickling requires to keep the object
        rng = read_psf_range(None, state)
        self.org = rng.org
        self.chr = rng.chr
        self.start = rng.start
        self.end = rng.end

    def __eq__(l, r):
        if not isinstance(r, Range):
            return False
        return l.org == r.org and l.chr == r.chr and \
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
        return hash(self.org) + hash(self.chr) + hash(self.start) + hash(self.end)

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
        return Range(self.org, self.chr, self.start + start, self.end - end)

    def is_inverted(self):
        return self.end < self.start

    def check(self):
        """
        DEPRECATED, since declaring positions as unsigned ints
        Does some basic sanity-checking if any positions are negative.
        """
        if self.start < 0:
            return False
        if self.end < 0:
            return False
        return True

cpdef read_psf_range(org:str, cell: str):
    """Parse a Range in PSF format
    PSF format contains the chromosome identifier followed by a colon (:) and the start and end position separated by -.
    Optionally, the organism may be specified as well, at the first position and separated by a colon (:). This overrides `org`, if it is passed.
    If the end position is before the start position, the range is treated as inverted.
    `org` specifies the organism the returned range should have. If this range is just used for filtering, None may be passed.
    Examples: Chr1:1000-2000, Chr3:10000-50000
    """
    #TODO error handling in here
    #print(cell)
    cellarr = cell.split(':')
    if len(cellarr) < 2 or len(cellarr) > 3:
        raise ValueError(f"Invalid PSF Range string: {cell}")
    if len(cellarr) == 3: # if a ref name is specified in the cell, that overrides the argument
        org = cellarr[0]
        cellarr = cellarr[1:]

    posarr = cellarr[1].split('-')
    cdef:
        str chrom = cellarr[0] #util.chrom_to_int(cellarr[0])
        int start = int(posarr[0])
        int end = int(posarr[1])
    return Range(org, chrom, start, end)


cpdef read_psf_pos(org:str, cell: str):
    """Parse a point position in PSF format
    PSF format contains the chromosome identifier followed by a colon (:) and the position.
    Optionally, the organism may be specified as well, at the first position and separated by a colon (:). This overrides `org`, if it is passed.
    Examples: Chr1:1000, Chr3:10000
    """
    cellarr = cell.split(':')
    if len(cellarr) < 2 or len(cellarr) > 3:
        raise ValueError(f"Invalid PSF Range string: {cell}")
    if len(cellarr) == 3: # if a ref name is specified in the cell, that overrides the argument
        org = cellarr[0]
        cellarr = cellarr[1:]

    return Position(org, cellarr[0], int(cellarr[1]))
