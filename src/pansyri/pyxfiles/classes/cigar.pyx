#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import itertools
import logging
import re
import copy

from cpython cimport array
import array

logger = logging.getLogger(__name__)

## constants

#TODO try making these Cpp sets
cdef reffwd = set(['M', 'D', 'N', '=', 'X'])
cdef qryfwd = set(['M', 'I', 'S', '=', 'X'])
cdef cig_types = set(['M', '=', 'X', 'S', 'H', 'D', 'I', 'N'])
cdef cig_aln_types = set(['M', 'X', '='])
cdef cig_clips = set(['S', 'H', 'P', 'N']) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these

cdef c_reffwd = set([ord('M'), ord('D'), ord('N'), ord('='), ord('X')])
cdef c_reffwd_noclip = set([ord('M'), ord('D'), ord('='), ord('X')])
cdef c_qryfwd = set([ord('M'), ord('I'), ord('S'), ord('='), ord('X')])
cdef c_qryfwd_noclip = set([ord('M'), ord('I'), ord('='), ord('X')])
cdef c_cig_types = set([ord('M'), ord('='), ord('X'), ord('S'), ord('H'), ord('D'), ord('I'), ord('N')])
cdef c_cig_aln_types = set([ord('M'), ord('X'), ord('=')])
cdef c_cig_clips = set([ord('S'), ord('H'), ord('P'), ord('N')]) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these

cdef relen = r"(\d+)[=XIDMNSHP]"
cdef retype = r"\d+([=XIDMNSHP])"

cdef inttr = lambda t: int(t[1])
cdef chrtr = lambda t: ord(t[1])

# declared outside of Cigar to be accessible from python, might move back later
cpdef cigar_from_string(str cg):
    """
    Takes a cigar string as input and returns a Cigar tuple
    """
    lens = array.array('I', [])
    types = array.array('b', [])
    
    # calling .extend should automatically preallocate when called by cython
    lens.extend(map(inttr, re.finditer(relen, cg)))
    types.extend(map(chrtr, re.finditer(retype, cg)))
    return Cigar.__new__(Cigar, lens, types)

# maybe implement cigar_from_full_string?

# try with the two array-based strategy again, try out if it's faster

cdef class Cigar:
    cdef array.array lens # array storing the lengths of each cigar tuple
    cdef array.array types # array storing the type of each cigar tuple

    def __cinit__(self, lens, types):
        self.lens = lens
        self.types = types

    def get_len(self, bint ref=True):
        """
        Returns the number of bases covered by this Cigar in the reference or query genome.
        """
        return self.get_len_of_type(c_reffwd_noclip if ref else c_qryfwd_noclip)

    def get_len_of_type(self, typeset):
        cdef unsigned int buf = 0
        for i in range(len(self.lens)):
            if self.types[i] in typeset:
                buf += self.lens[i]
        return buf

    def get_identity(self, bint ref=True):
        """
        Returns the fraction of covered bases (of the reference/query) that are an exact match ('=').
        """
        return self.get_len_of_type(set(ord('=')))/self.get_len(ref=ref)

    def __len__(self):
        return sum(self.lens)

    def __eq__(self, other):
        if isinstance(other, Cigar):
            return self.lens == other.lens and self.types == other.types
        else:
            return False

    def __ne__(self, other):
        if isinstance(other, Cigar):
            return not(self.lens == other.lens and self.types == other.types)
        else:
            return True

    def pad(self, unsigned int left, unsigned int right, clip='S'):
        """Small function that adds padding to one or both sides of this `Cigar`.
        Mutates self!

        :param left, right: How much padding to add to each side.
        :type left, right: `int`
        :param clip: which kind of padding to use, can be 'S' or 'H'. 'S' by default.
        :return: `None`
        """
        if left > 0:
            self.lens.insert(0, left)
            self.types.insert(0, ord(clip))
        if right > 0:
            self.lens.append(right)
            self.types.append(0, ord(clip))

    def unpad(self):
        """Removes padding from the borders of this `Cigar`.
        Mutates self!
        :return: `None`
        """
        cdef size_t start = 0
        for tup in self.tups:
            if tup.t not in c_cig_clips: # wtf is up here? why is this converted to a dict, but this is not done in get_len_of_type?
                # must be something with the context? no idea
                break
            start += 1

        self.tups.erase(self.tups.begin(), self.tups.begin() + start)

        cdef size_t i_end = 0
        for tup in self.tups[::-1]:
            if tup.t not in cig_clips:
                break
            i_end += 1

        self.tups.erase(self.tups.end()-i_end, self.tups.end())

    def reverse(self):
        """
        Returns a Cigar representing a reversed version of the alignment represented by this Cigar.
        """
        newlens = copy.deepcopy(self.lens)
        newtypes = copy.deepcopy(self.types)
        newlens.reverse()
        newtypes.reverse()
        return Cigar(newlens, newtypes)

    # TODO maybe benchmark bints vs separate char's or argstruct or separate methods
    cpdef get_removed(self, unsigned int n, bint ref=True, bint start=True, bint only_pos=False):
        """
        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
        """
        if len(self.lens) == 0:
            logger.error("Trying to remove from an empty Cigar!")
            raise ValueError("empty Cigar!")

        if n == 0: # shortcut for a common path
            if only_pos:
                return 0
            else:
                return (0, self)
        
        cdef:
            int ind = 0 # position currently being evaluated for skipping
            unsigned int skip = 0 # bases skipped in the other sequence
            # two sets containing the CIGAR codes incrementing one or the other strand
            fwd = c_reffwd if ref else c_qryfwd 
            altfwd = c_qryfwd if ref else c_reffwd
            int rem = n # tally how much is still left to remove
            curt = self.types[ind] if start else self.types[-1]
            curn = self.lens[ind] if start else self.lens[-1]

        # loop and remove regions as long as the skip is more than one region
        while rem > 0 and ind < len(self.types):
            curt = self.types[ind] if start else self.types[-ind-1]
            curn = self.lens[ind] if start else self.lens[-ind-1]
            # increment appropriate counters depending on which strand this cgi forwards
            if curt in altfwd:
                skip += curn
            if curt in fwd:
                rem -= curn
            ind += 1

        if rem > 0:
            logger.error(f"tried to remove more than CIGAR length Params: n: {n}, start: {start}, ref: {ref}, Cigar length: {len(self)}, terminated at index {ind}")
            raise ValueError("tried to remove more than CIGAR length")

        if curt in altfwd: # remove overadded value
            skip += rem

        if only_pos:
            return skip
        
        newlens = self.lens[ind:] if start else self.lens[:-ind]
        newtypes = self.types[ind:] if start else self.types[:-ind]

        # if there is a remainder, add it to the front
        if rem < 0:
            if start:
                newlens.insert(0, -rem)
                newtypes.insert(0, curt)
            else:
                newlens.append(-rem)
                newtypes.append(curt)

        return (skip, Cigar(newlens, newtypes))


    def __repr__(self):
        return f"Cigar({self.to_string()})"

    def to_string(self):
        return ''.join([str(tup.n) + chr(tup.t) for tup in self.tups])

    def to_full_string(self):
        return ''.join([chr(tup.t)*tup.n for tup in self.tups])


#TODO rewrite with copying for cython
#    def clean(self):
#        """
#        Misc function to remove empty annotations and combine neighbouring annotations of equal type from a Cigar.
#        Mutates self, but the represented alignment stays the same.
#        """
#        i = 1
#        while i < len(self.pairs):
#            if self.pairs[i-1][0] == 0:
#                del self.pairs[i-1]
#                continue
#            if self.pairs[i-1][1] == self.pairs[i][1]:
#                self.pairs[i-1][0] += self.pairs[i][0]
#                del self.pairs[i]
#            else:
#                i += 1
