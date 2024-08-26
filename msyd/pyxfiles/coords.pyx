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
#TODO use proper cython types, e.g. char for haplo

# decorator to auto-implement __gt__ etc. from __lt__ and __eq__
@functools.total_ordering # not sure how performant, TO/DO replace later?
cdef class Position:
    def __init__(self, org:str, chr:int, haplo:str, pos: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.pos = pos

    def __repr__(self):
        #return f"Position({self.org}, {self.chr}, {self.haplo}, {self.pos})"
        return f"Position({self.org}, {self.chr}, {self.pos})"

    def to_psf(self):
        """Transform this `Position` into population synteny file format
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
cdef class Range:
    def __init__(self, org:str, chr:str, haplo, start: int, end: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.start = start # inclusive
        self.end = end # inclusive

    def __repr__(self):
        #return f"Range({self.org}, {self.chr}, {self.haplo}, {self.start}, {self.end})"
        return f"Range({self.org}, {self.chr}, {self.start}, {self.end})"

    def to_psf(self):
        """Transform this `Range` into the form specified by PSF
        """
        #return f"{self.chr}:{self.haplo}:{self.start}-{self.end}"
        return f"{self.chr}:{self.start}-{self.end}"

    def to_psf_org(self):
        """Transform this `Range` into the form specified by PSF, with the sample name being prepended as specified for the realigned reference haplotype
        """
        #return f"{self.org}:{self.chr}:{self.haplo}:{self.start}-{self.end}"
        return f"{self.org}:{self.chr}:{self.start}-{self.end}"


    def read_psf(org:str, cell: str):
        """Parse a Range in PSF format
        PSF format is :-separated and must contain the chromosome first, then followed by a haplotype (optional) and a start and end position separated by -.
        If the end position is before the start position, the range is treated as inverted.
        `org` specifies the organism the returned range should have. If this range is just used for filtering, None may be passed.
        Examples: Chr1:mat:1000-2000, Chr3:10000-50000
        """
        #TODO error handling in here
        #print(cell)
        cellarr = cell.split(':')
        if len(cellarr) < 2 or len(cellarr) > 4:
            raise ValueError(f"Invalid PSF Range string: {cell}")
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
