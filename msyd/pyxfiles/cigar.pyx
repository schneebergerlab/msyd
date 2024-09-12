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
from libcpp.unordered_set cimport unordered_set

import msyd.util as util


logger = util.CustomFormatter.getlogger(__name__)

## constants

#TODO try making these Cpp sets
cdef:
    reffwd = set(['M', 'D', 'N', '=', 'X'])
    qryfwd = set(['M', 'I', 'S', '=', 'X'])
    cig_types = set(['M', '=', 'X', 'S', 'H', 'D', 'I', 'N'])
    cig_aln_types = set(['M', 'X', '='])
    cig_clips = set(['S', 'H', 'P', 'N']) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these

    unordered_set[char] c_reffwd = unordered_set[char]([ord('M'), ord('D'), ord('N'), ord('='), ord('X')])
    unordered_set[char] c_reffwd_noclip = unordered_set[char]([ord('M'), ord('D'), ord('='), ord('X')])
    unordered_set[char] c_qryfwd = unordered_set[char]([ord('M'), ord('I'), ord('S'), ord('='), ord('X')])
    unordered_set[char] c_qryfwd_noclip = unordered_set[char]([ord('M'), ord('I'), ord('='), ord('X')])
    unordered_set[char] c_cig_types = unordered_set[char]([ord('M'), ord('='), ord('X'), ord('S'), ord('H'), ord('D'), ord('I'), ord('N')])
    unordered_set[char] c_cig_aln_types = unordered_set[char]([ord('M'), ord('X'), ord('=')])
    unordered_set[char] c_cig_clips = unordered_set[char]([ord('S'), ord('H'), ord('P'), ord('N')]) # N is not clipping, but is ignored anyway. Really, it shouldnord('t even occur in alignments like these

    bam_code_map = [ord('M'), ord('I'), ord('D'), ord('N'), ord('S'), ord('H'), ord('P'), ord('='), ord('X')]

    retup = r"(\d+)([=XIDMNSHP])"

# declared outside of Cigar to be accessible from python, might move back later
cpdef Cigar cigar_from_string(str cg):
    """
    Takes a cigar string as input and returns a Cigar object
    """
    return Cigar.__new__(Cigar, tups=cigt_from_string(cg))

cdef vector[Cigt] cigt_from_string(str cg):
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
    
    # free up unnecessarily reserved memory in case we were too pessimistic
    tups.shrink_to_fit()

    return tups

# maybe implement cigar_from_full_string?
cpdef cigar_from_bam(bam):    
    """
    Takes a List of Cigar tuples with BAM codes as input, returns as a Cigar struct.
    """
    cdef vector[Cigt] tups = vector[Cigt]()
    # preallocate assuming on average each tuple has two digit length
    # maybe try being more optimistic and assuming three-digit length
    tups.reserve(len(bam))

    for tup in bam:
        assert(0 < tup[1] and tup[1] < 9)
        tups.push_back(Cigt(tup[0], bam_code_map[tup[1]]))
 
    return Cigar.__new__(Cigar, tups=tups)

# small struct to contain the length and type of a cigar tuple
cdef packed struct Cigt:
#cdef struct Cigt: # slower
    unsigned int n
    char t



# got it working to not work with two arrays
# pretty sure this is faster, might try exact benchmark though


cdef class Cigar:
    #cdef array.array lens # array storing the lengths of each cigar tuple
    #cdef array.array types # array storing the type of each cigar tuple
    cdef vector[Cigt] tups

    def __cinit__(self, tups=None):
        """
        tups is optional to support pickling; always set it otherwise!
        """
        if tups is not None:
            self.tups = tups
        else:
            self.tups = vector[Cigt]()

    def get_len(self, bint ref=True):
        """
        Returns the number of bases covered by this Cigar in the reference or query genome.
        """
        return self.get_len_of_type(c_reffwd_noclip if ref else c_qryfwd_noclip)

    cdef get_len_of_type(self, unordered_set[char] typeset):
        cdef unsigned int buf = 0
        for tup in self.tups:
            if typeset.count(tup.t): # contains method still not supported until C++20
                buf += tup.n
        return buf

    def get_identity(self):#, bint ref=True):
        """
        Returns the fraction of covered bases (of the reference/query) that are an exact match ('=').
        """
        return self.get_len_of_type(unordered_set[char]({ord('=')}))/len(self)

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

    def __repr__(self):
        return f"Cigar({self.to_string()})"

    def to_string(self):
        return ''.join([str(tup.n) + chr(tup.t) for tup in self.tups])

    def to_full_string(self):
        return ''.join([chr(tup.t)*tup.n for tup in self.tups])

    # support pickling, for use with multiprocessing
    def __getstate__(self):
        return self.to_string()

    def __setstate__(self, state):
        self.tups = cigt_from_string(state)

    def is_empty(self):
        return self.tups.empty()

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
            if not c_cig_clips.count(tup.t): # wtf is up here? why is this converted to a dict, but this is not done in get_len_of_type?
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
        cdef vector[Cigt] newtups = vector[Cigt]()
        newtups.reserve(self.tups.size())
        for i in range(self.tups.size()):
            newtups.push_back(self.tups[self.tups.size() -i -1])

        return Cigar(newtups)

    cpdef trim(self, unsigned int s, unsigned int e, bint ref=True, bint only_pos=False):
        """
        Trims an alignment by removing `s` bases from the start and `e` from the end.
        If `ref` is set to `True`, the removed bases are counted on the reference sequence, otherwise on the alternative.
        Returns a tuple with the first value representing the number of bases removed from the start in the alternative (resp. reference if `ref` is `True`).
        The second value is the same, but from the end of the sequence.
        If `only_pos` is set, will not change the Cigar, but only return the number of bases that would be removed.
        Internally calls get_removed().
        """
        if only_pos:
            return (self.get_removed(s, ref=ref, only_pos=only_pos), self.get_removed(e, ref=ref, only_pos=only_pos))
        sdrop, tmp = self.get_removed(s, ref=ref)
        edrop, tmp = tmp.get_removed(e, ref=ref)
        return (sdrop, edrop, tmp)

    cpdef trim_matching(self):
        """
        Trims a CIGAR string until both ends start with a matching (=) position.
        :returns: The number of bases deleted in the query/ref and a new CIGAR guaranteed to start and end with =.
        """
        # trim from start until = found
        cdef:
            int start = 0
            int qstart = 0
            int rstart = 0
        cdef Cigt cur = self.vector[start]
        while cur.t != ord('='):
            if c_reffwd.count(cur.t): # skip on r
                rstart += cur.n
            if c_qryfwd.count(cur.t): # skip on q
                qstart += cur.n
            # go forward in the loop
            start += 1
            cur = self.vector[start]

        # trim from end until = found
        cdef:
            int end = self.vector.end()
            int qend = 0
            int rend = 0
        cur = self.vector[end]
        while cur.t != ord('='):
            if c_reffwd.count(cur.t): # skip on r
                rend += cur.n
            if c_qryfwd.count(cur.t): # skip on q
                qend += cur.n
            # go forward in the loop
            end -= 1
            cur = self.vector[end]

        # construct new cigar using indexes
        cdef vector[Cigt] newcg = vector[Cigt]()
        newcg.reserve(end - start)
        # add Cigts to newcg
        while start <= end:
            newcg.push_back(self.vector[start])
            start +=1

        return qstart, qend, rstart, rend, Cigar(newcg)



    # TODO maybe benchmark bints vs separate char's or argstruct or separate methods
    cpdef get_removed(self, unsigned int n, bint ref=True, bint start=True, bint only_pos=False): #nogil
        """
        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
        """

        if n == 0: # shortcut for a common path
            if only_pos:
                return 0
            else:
                return (0, self)
        # deleting from an empty record is fine, so long as 0 is deleted
        if self.tups.empty():
            logger.error("Trying to remove from an empty Cigar!")
            raise ValueError("empty Cigar!")
        
        cdef:
            unsigned int ind = 0 # position currently being evaluated for skipping
            unsigned int skip = 0 # bases skipped in the other sequence
            # two sets containing the CIGAR codes incrementing one or the other strand
            unordered_set[char] fwd = c_reffwd if ref else c_qryfwd 
            unordered_set[char] altfwd = c_qryfwd if ref else c_reffwd
            int rem = n # tally how much is still left to remove
            Cigt cur = self.tups[ind] if start else self.tups[self.tups.size()-1]

        # loop and remove regions as long as the skip is more than one region
        # >= to greedily remove I/D at the edges
        while rem >= 0 and ind < self.tups.size():
            cur = self.tups[ind] if start else self.tups[self.tups.size()-ind -1]
            # increment appropriate counters depending on which strand this cgi forwards
            if altfwd.count(cur.t):
                skip += cur.n
            if fwd.count(cur.t):
                rem -= cur.n
            ind += 1

        if rem > 0:
            logger.error(f"tried to remove more than CIGAR length Params: n: {n}, start: {start}, ref: {ref}, Cigar len on ref/alt: {self.get_len(ref=ref)}, terminated at index {ind}")
            raise ValueError("tried to remove more than CIGAR length")

        if altfwd.count(cur.t): # remove overadded value
            skip += rem

        if only_pos:
            return skip
        # TODO try changing backing vector and only storing index
        
        cdef vector[Cigt] newtups = vector[Cigt]()
        newtups.reserve(self.tups.size() - ind + 1)
        if start:
            # if there is a remainder, add it to the front
            if rem < 0:
                newtups.push_back(Cigt(-rem, cur.t))
            # TODO check if this is cythonized efficiently 
            for i in range(ind, self.tups.size()):
                newtups.push_back(self.tups[i])
        else:
            for i in range(self.tups.size()-ind):
                newtups.push_back(self.tups[i])
            # if there is a remainder, add it to the back
            if rem < 0:
                newtups.push_back(Cigt(-rem, cur.t))

        return (skip, Cigar(newtups))



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
