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
        return sum([i[0] for i in self.pairs if i[1] in s])

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
        # TODO proper error handling when removing more than possible or if input null or sth 
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
        

        # loop and remove regions as long as the skip is more than one region
        try:
            while ind < len(self.pairs):
                cgi = self.pairs[ind] if start else self.pairs[-ind -1]
                #print(ind, cgi, n, skip)
                if cgi[1] not in fwd:
                    ind += 1
                    if cgi[1] in altfwd:
                        skip += cgi[0]
                    continue
                
                # this region counts towards n, determine if it can be removed or is too large
                if n >= cgi[0]:
                    n -= cgi[0]
                    ind += 1
                    if cgi[1] in altfwd:
                        skip += cgi[0]
                else:
                    break
        except IndexError:
            print("ERROR: not removing more than sequence length, returning None")
            raise ValueError
            return None

        # remaining skip must be < 1 in a fwd-region
        # construct the new return cigar, without altering self
        cgi = self.pairs[ind] if start else self.pairs[-ind -1]
        if cgi[1] in altfwd:
            skip += n

        if start:
            return (skip, Cigar([[cgi[0]-n, cgi[1]]] + self.pairs[ind+1:]))
        else:
            return (skip, Cigar(self.pairs[:-ind-1] + [[cgi[0]-n, cgi[1]]]))

    def __repr__(self):
        return f"Cigar({self.pairs})"

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


    def impute(l, r):
        """
        This function combines two CIGAR strings of two queries on one reference and tries to impute the CIGAR string for a 1:1 alignment of one query to the other.
        By necessity, this is a crude approximation and can in no way replace a full alignment.
        The algorithm assumes the alignments to not start shifted from each other.
        :param: Two Cigars.
        :return: A Cigar.
        """
        imputed_pairs = []
        if len(l.pairs) == 0 or len(r.pairs) ==0:
            raise ValueError("Empty Cigar input into impute!")

        # prepare the variables for iteration
        liter = iter(l.pairs)
        riter = iter(r.pairs)
        lpr = next(liter)
        rpr = next(riter)
        lprog = 0
        rprog = 0

        while(True):
            try:
                # check if a region has been fully consumed
                if lpr[0]-lprog <= 0 or lpr[1] in cig_clips:
                    lprog = 0
                    lpr = next(liter)
                if rpr[0]-rprog <= 0 or rpr[1] in cig_clips:
                    rprog = 0
                    rpr = next(riter)

                # compute the max possible step
                step = min(lpr[0]-lprog, rpr[0]-rprog)
                # compute the type and which strand to step
                t, stepl, stepr = Cigar.impute_type(lpr[1], rpr[1])
                if t is not None:
                    imputed_pairs.append([step, t])
                # do the step #TO/DO make this branchless with cdefs
                if stepr:
                    rprog += step
                if stepl:
                    lprog += step

            except StopIteration:
                break

        # append clippings at end for the sequence that hasn't ran out
        #for p in liter:
        #    imputed_pairs.append([p[0], 'S'])
        #for p in riter:
        #    imputed_pairs.append([p[0], 'S'])

        # return
        cg = Cigar(imputed_pairs)
        cg.clean()
        return cg

    def impute_type(ltp, rtp):
        """
        Helper function for impute().
        Given two CIGAR types, outputs the one that is imputed
        :param: two valid CIGAR chars (M=IDSHX)
        :return: the CIGAR char imputed for the inputs as well as two bools specifying whether to step ahead on the left/right Cigar.
        """
        if ltp in cig_aln_types:
            if rtp == 'D':
                return 'D', True, True

            if rtp == 'I':
                return 'I', True, True

            if ltp == 'X' or rtp == 'X':
                return 'X', True, True
            else:
                if ltp == '=' and rtp == '=':
                    return '=', True, True
                else: # at least one must be an M, none is an X
                    return 'M', True, True

        if ltp == 'I':
            if rtp == 'I': # TO/DO be even more conservative and resolve this to D followed by I?
                return 'M', True, True
            else:
                return 'D', True, False # right has deletion relative to left
            
        if ltp == 'D':
            if rtp == 'D':
                return None, True, True # None to skip this region
            elif rtp in cig_aln_types:
                return 'I', True, True
            elif rtp == 'I':
                return 'I', True, True


        


if __name__ == "__main__":
    import sys

    for x in sys.argv[1:]:
        cg = Cigar.from_string(x)
        print(cg)
        print(cg.get_removed_faster(10), cg.get_removed(10), sep='\n')
        print(cg.get_removed_faster(10, ref=False), cg.get_removed(10, ref=False), sep='\n')
        print(cg.get_removed_faster(10, start=False), cg.get_removed(10, start=False), sep='\n')
        print(cg.get_removed_faster(10, start=False, ref=False), cg.get_removed(10, start=False, ref=False), sep='\n')

        
    

    sys.exit()
    import ingest
    
    print(Cigar.from_full_string(sys.argv[1]).to_string())
    sys.exit()

    dfr = ingest.readSAMBAM(sys.argv[1])
    dfq = ingest.readSAMBAM(sys.argv[2])
    dfc = ingest.readSAMBAM(sys.argv[3])

    rowr = dfr.loc[0]
    rowq = dfq.loc[1]
    rowc = dfc.loc[1]

    print(rowr)
    print(rowq)
    print(rowc)


    # find the overlap
    start = max(rowr['astart'], rowq['astart'], rowc['astart'])
    end = min(rowr['aend'], rowq['aend'], rowc['aend'])
    print(start, end)

    strskpr, endskpr = start - rowr['astart'], rowr['aend'] - end
    strskpq, endskpq = start - rowq['astart'], rowq['aend'] - end
    strskpc, endskpc = start - rowc['astart'], rowc['aend'] - end

    # QUERY sequences
    sqr = dfr.loc[0, 'seq'][strskpr:-endskpr-1]
    sqq = dfq.loc[1, 'seq'][strskpq:-endskpq-1]
    sqc = dfc.loc[1, 'seq'][strskpc:-endskpc-1]

    cgr = Cigar.from_string(dfr.loc[0, 'cg']).get_removed(strskpr)[1].get_removed(endskpr, start=False)[1]
    cgq = Cigar.from_string(dfq.loc[1, 'cg']).get_removed(strskpq)[1].get_removed(endskpq, start=False)[1]
    cgc = Cigar.from_string(dfc.loc[1, 'cg']).get_removed(strskpc)[1].get_removed(endskpc, start=False)[1]
    
    assert(cgr.get_len() == cgq.get_len() == cgc.get_len())

    cgi = Cigar.impute(cgr, cgq)

    # prettyprint
    cgcstr = cgc.to_full_string()
    cgistr = cgi.to_full_string()
    #for n in range(79, len(cgcstr), 79):
        #if all([cgistr[i]==cgcstr[i] and cgcstr[i]=='=' for i in range(n-79,n)]):
        #    continue
    #    print(">", cgcstr[n-79:n])
    #    print("<", cgistr[n-79:n])


