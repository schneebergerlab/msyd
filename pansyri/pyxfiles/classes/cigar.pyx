#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import re
import copy
import itertools

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

    def __eq__(l, r):
        return l.pairs == r.pairs

    def get_identity(self):
        """
        Returns the fraction of covered bases (of the reference) that are an exact match ('=').
        """
        return sum([int(i[0]) for i in self.pairs if i[1] == '='])/self.get_len()

    def get_removed_legacy(self, n, ref=True, start=True):
        """
        If ref=True, removes from the 'start'/end of the QUERY strand until 'n' bases from the REFERENCE strand have been removed, if ref=False vice versa.
        :return: The number of bases deleted in the query/ref and a CIGAR with these bases removed.
        """
        cg = copy.deepcopy(self)

        ind = 0 if start else -1
        skip = 0
        fwd = reffwd if ref else qryfwd
        altfwd = qryfwd if ref else reffwd

        while n > 0:
            #TODO speed this up by first determining the index to subset,
            # then copy the subset and finally adjust the border element
            # TODO implement this
            cgi = cg.pairs[ind]
            #print(n, start, cgi)
            if cgi[1] in altfwd:
                skip += cgi[0]

            if cgi[1] not in fwd:
                # if ever removing the copy, refactor!
                del cg.pairs[ind]
                continue

            # cgi must be in fwd
            n -= cgi[0]
            if n >= 0:
                del cg.pairs[ind]
            else:
                cgi[0] = -n # n == n- cg[ind][0] implies -n == cg[ind][0] -n
                if cgi[1] in altfwd: # subtract the overcounting
                    skip += n

        return (skip, cg)


    def get_removed(self, n, ref=True, start=True):
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
        rem = n # tally how much still left to remove
        

        # loop and remove regions as long as the skip is more than one region
        try:
            while ind < len(self.pairs):
                cgi = self.pairs[ind] if start else self.pairs[-ind -1]
                if cgi[1] not in fwd:
                    ind += 1
                    if cgi[1] in altfwd:
                        skip += cgi[0]
                # this region counts towards n, determine if it can be removed or is too large
                elif rem >= cgi[0]:
                    rem -= cgi[0]
                    ind += 1
                    if cgi[1] in altfwd:
                        skip += cgi[0]
                else:
                    break

            cgi = self.pairs[ind] if start else self.pairs[-ind -1]
            if cgi[1] in altfwd:
                skip += rem

            if start:
                return (skip, Cigar([[cgi[0]-rem, cgi[1]]] + self.pairs[ind+1:]))
            else:
                return (skip, Cigar(self.pairs[:-ind-1] + [[cgi[0]-rem, cgi[1]]]))

        except IndexError:
            print(f"ERROR: tried to remove more than sequence length in get_removed of {n} with start {start} on ref {ref} on Cigar with length {len(self)} at index {ind}")
            raise ValueError("invalid skip")

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