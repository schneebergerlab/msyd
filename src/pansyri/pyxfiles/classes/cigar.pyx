#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import itertools
import logging
import re

from cpython cimport array
import array

from libcpp.vector cimport vector


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
cdef c_cig_clips = set([ord('S'), ord('H'), ord('P'), ord('N')]) # N is not clipping, but is ignored anyway. Really, it shouldnord('t even occur in alignments like these

cdef retup = r"(\d+)([=XIDMNSHP])"

# declared outside of Cigar to be accessible from python, might move back later
cpdef cigar_from_string(str cg):
    """
    Takes a cigar string as input and returns a Cigar tuple
    """
    cdef vector[Cigt] tups = vector[Cigt]()
    # preallocate assuming on average each tuple has two digit length
    # maybe try being more optimistic and assuming three-digit length
    tups.reserve(int(len(cg)/3))

    for match in re.findall(retup, cg):
        if match[1] not in cig_types:
            logger.error("Tried to construct a Cigar object with invalid type")
            raise ValueError("Not a CIGAR type!")
        tups.push_back(Cigt(int(match[0]), ord(match[1])))
 
    return Cigar.__new__(Cigar, tups)

# maybe implement cigar_from_full_string?
    
# small struct to contain the length and type of a cigar tuple
cdef packed struct Cigt:
#cdef struct Cigt: # slower
    unsigned short n
    char t



# got it working to not work with two arrays
# pretty sure this is faster, might try exact benchmark though


cdef class Cigar:
    #cdef array.array lens # array storing the lengths of each cigar tuple
    #cdef array.array types # array storing the type of each cigar tuple
    cdef vector[Cigt] tups

    def __cinit__(self, tup):
        self.tups = tup

    def get_len(self, bint ref=True):
        """
        Returns the number of bases covered by this Cigar in the reference or query genome.
        """
        return self.get_len_of_type(c_reffwd_noclip if ref else c_qryfwd_noclip)

    def get_len_of_type(self, typeset):
        cdef unsigned int buf = 0
        for tup in self.tups:
            if tup.t in typeset:
                buf += tup.n
        return buf

    def get_identity(self, bint ref=True):
        """
        Returns the fraction of covered bases (of the reference/query) that are an exact match ('=').
        """
        return self.get_len_of_type(set(ord('=')))/self.get_len(ref=ref)

    def __len__(self):
        cdef unsigned int buf = 0
        for tup in self.tups:
            buf += tup.n
        return buf

    # internal cdef'd method to avoid exposing tups directly to python
    cdef equals(self, Cigar other):
        if self.tups.size() != other.tups.size():
            return False
        for i in range(self.tups.size()):
            selft = self.tups.at(i)
            othert = other.tups.at(i)
            # apparently, structs cannot have associated methods or overloaded operators
            if selft.t != othert.t or selft.n != othert.n:
                return False
        return True

    def __eq__(self, other):
        if isinstance(other, Cigar):
            return self.equals(other)
        else:
            return False

    def __ne__(self, other):
        if isinstance(other, Cigar):
            return not self.equals(other)
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
            p = Cigt(left, ord(clip))
            self.tups.insert(self.tups.begin(), p)
        if right > 0:
            self.tups.push_back(Cigt(right, ord(clip)))

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


    # TODO maybe benchmark bints vs separate char's or argstruct or separate methods
    cpdef get_removed(self, unsigned int n, bint ref=True, bint start=True, bint only_pos=False):
        """
        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
        """
        if self.tups.empty():
            logger.error("Trying to remove from an empty Cigar!")
            raise ValueError("empty Cigar!")

        if n == 0: # shortcut for a common path
            if only_pos:
                return 0
            else:
                return (0, self)
        
        cdef:
            unsigned int ind = 0 # position currently being evaluated for skipping
            unsigned int skip = 0 # bases skipped in the other sequence
            # two sets containing the CIGAR codes incrementing one or the other strand
            fwd = c_reffwd if ref else c_qryfwd 
            altfwd = c_qryfwd if ref else c_reffwd
            int rem = n # tally how much is still left to remove
            Cigt cur = self.tups[ind] if start else self.tups[self.tups.size()-1]

        # loop and remove regions as long as the skip is more than one region
        while rem > 0 and ind < self.tups.size():
            cur = self.tups[ind] if start else self.tups[self.tups.size()-ind -1]
            # increment appropriate counters depending on which strand this cgi forwards
            if cur.t in altfwd:
                skip += cur.n
            if cur.t in fwd:
                rem -= cur.n
            ind += 1

        if rem > 0:
            logger.error(f"tried to remove more than CIGAR length Params: n: {n}, start: {start}, ref: {ref}, Cigar length: {len(self)}, terminated at index {ind}")
            raise ValueError("tried to remove more than CIGAR length")

        if cur.t in altfwd: # remove overadded value
            skip += rem

        if only_pos:
            return skip
        
        cdef vector[Cigt] newtups = vector[Cigt]()
        newtups.reserve(self.tups.size() - ind + 1)
        if start:
            # if there is a remainder, add it to the front
            if rem < 0:
                newtups.push_back(Cigt(-rem, cur.t))
            for i in range(ind, self.tups.size()):
                newtups.push_back(self.tups[i])
        else:
            for i in range(self.tups.size()-ind):
                newtups.push_back(self.tups[i])
            # if there is a remainder, add it to the back
            if rem < 0:
                newtups.push_back(Cigt(-rem, cur.t))

        return (skip, Cigar(newtups))


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
