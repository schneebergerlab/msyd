#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import re
import itertools
import logging

logger = logging.getLogger(__name__)

## constants
reffwd = set(['M', 'D', 'N', '=', 'X'])
qryfwd = set(['M', 'I', 'S', '=', 'X'])
cig_types = set(['M', '=', 'X', 'S', 'H', 'D', 'I', 'N'])
cig_aln_types = set(['M', 'X', '='])
cig_clips = set(['S', 'H', 'P', 'N']) # N is not clipping, but is ignored anyway. Really, it shouldn't even occur in alignments like these
inttr = lambda x: [int(x[0]), x[1]]

class Cigar:
    """
    A CIGAR string represented as a list of lists.
    """

    def __init__(self, pairs):
        self.pairs = pairs

    ## copied from myUsefulFunctions.py
    def from_string(cg):
        """
        Takes a cigar string as input and returns a Cigar tuple
        """
        #TODO more error handling
        
        for i in "MIDNSHPX=":
            cg = cg.replace(i, ';'+i+',')
        return Cigar([inttr(i.split(';')) for i in cg.split(',')[:-1]])

    def get_len(self, ref=True):
        """
        Returns the number of bases covered by this Cigar in the reference or query genome.
        """
        s = reffwd if ref else qryfwd
        return sum([i[0] for i in self.pairs if i[1] in s and not i[1] in cig_clips])

    def __len__(self):
        return sum([i[0] for i in self.pairs])

    def __eq__(self, other):
        if isinstance(other, Cigar):
            return self.pairs == other.pairs
        else:
            return False

    def __ne__(self, other):
        if isinstance(other, Cigar):
            return self.pairs != other.pairs
        else:
            return True

    def get_identity(self):
        """
        Returns the fraction of covered bases (of the reference) that are an exact match ('=').
        """
        return sum([int(i[0]) for i in self.pairs if i[1] == '='])/self.get_len()

    def pad(self, left: int, right: int, clip='S'):
        """Small function that adds padding to one or both sides of this `Cigar`.
        Mutates self!

        :param left, right: How much padding to add to each side.
        :type left, right: `int`
        :param clip: which kind of padding to use, can be 'S' or 'H'. 'S' by default.
        :return: `None`
        """
        if left < 0 or right < 0:
            raise ValueError("pad supplied with negative length!")

        if left > 0:
            self.pairs = [[left, clip]] + self.pairs
        if right > 0:
            self.pairs = self.pairs + [[right, clip]]

    def unpad(self):
        """Removes padding from the borders of this `Cigar`.
        Mutates self!
        :return: `None`
        """
        i_start = 0
        for x in self.pairs:
            if x[1] not in cig_clips:
                break
            i_start += 1
        i_end = len(self.pairs)
        for x in self.pairs[::-1]:
            if x[1] not in cig_clips:
                break
            i_end -= 1

        self.pairs = self.pairs[i_start:i_end]

    def unpad_all(self):
        """Removes all padding from this `Cigar`, even from the middle of the alignment (however it got there...)
        Mutates self!
        :return: `None`
        """
        for pos, pair in enumerate(self.pairs):
            if pair[1] in cig_clips:
                del self.pairs[pos]



    def get_removed(self, n, ref=True, start=True, only_pos=False):
        """
        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
        """
        if not self.pairs:
            raise ValueError("empty Cigar!")

        ind = 0 # position currently being evaluated for skipping
        skip = 0 # bases skipped in the other sequence
        # two sets containing the CIGAR codes incrementing one or the other strand
        fwd = reffwd if ref else qryfwd 
        altfwd = qryfwd if ref else reffwd
        rem = n # tally how much is still left to remove
        

        # loop and remove regions as long as the skip is more than one region
        while ind < len(self.pairs):
            cgi = self.pairs[ind] if start else self.pairs[-ind -1]
            # increment appropriate counters depending on which strand this cgi forwards
            if cgi[1] in altfwd:
                skip += cgi[0]
            if cgi[1] in fwd:
                rem -= cgi[0]
            # abort routine
            if rem <= 0:
                if cgi[1] in altfwd: # add remainder if necessary
                    skip += rem
                break
            ind += 1

        if rem > 0:
            logger.error(f"tried to remove more than CIGAR length Params: n: {n}, start: {start}, ref: {ref}, Cigar length: {len(self)}, terminated at index {ind}")
            raise ValueError("tried to remove more than CIGAR length")

        if only_pos:
            return skip

        newtuplist = [[-rem, cgi[1]]] if rem < 0 else []
        if start:
            return (skip, Cigar(newtuplist + self.pairs[ind+1:]))
        else:
            return (skip, Cigar(self.pairs[:-ind-1] + newtuplist))


    def __repr__(self):
        return f"Cigar({self.to_string()})"

    def to_string(self):
        return ''.join([str(p[0]) + p[1] for p in self.pairs])

    def to_full_string(self):
        return ''.join([p[0]*p[1] for p in self.pairs])

    def from_full_string(string):
        pairs = []
        mode = string[0]
        count = 0
        for ch in string:
            if ch == mode:
                count += 1
            else:
                pairs.append([count, mode])
                count = 1
                mode = ch

        pairs.append([count, mode])

        return Cigar(pairs)



    def clean(self):
        """
        Misc function to remove empty annotations and combine neighbouring annotations of equal type from a Cigar.
        Mutates self, but the represented alignment stays the same.
        """
        i = 1
        while i < len(self.pairs):
            if self.pairs[i-1][0] == 0:
                del self.pairs[i-1]
                continue
            if self.pairs[i-1][1] == self.pairs[i][1]:
                self.pairs[i-1][0] += self.pairs[i][0]
                del self.pairs[i]
            else:
                i += 1

#import copy
#
#    def get_removed_legacy(self, n, ref=True, start=True):
#        """
#        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
#        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
#        """
#        cg = copy.deepcopy(self)
#
#        ind = 0 if start else -1
#        skip = 0
#        fwd = reffwd if ref else qryfwd
#        altfwd = qryfwd if ref else reffwd
#
#        while n > 0:
#            #TODO speed this up by first determining the index to subset,
#            # then copy the subset and finally adjust the border element
#            # TODO implement this
#            cgi = cg.pairs[ind]
#            #logger.info(f"{n}, {start}, {cgi}")
#            if cgi[1] in altfwd:
#                skip += cgi[0]
#
#            if cgi[1] not in fwd:
#                # if ever removing the copy, refactor!
#                del cg.pairs[ind]
#                continue
#
#            # cgi must be in fwd
#            n -= cgi[0]
#            if n >= 0:
#                del cg.pairs[ind]
#            else:
#                cgi[0] = -n # n == n- cg[ind][0] implies -n == cg[ind][0] -n
#                if cgi[1] in altfwd: # subtract the overcounting
#                    skip += n
#
#        return (skip, cg)
