# -*- coding: utf-8 -*-
#!/usr/bin/python3
# distutils: language = c++
# cython: language_level = 3

"""
This file will find pansyntenic regions, either to use in multiple structural variant calling or for separate analysis.
"""
import functools
import pandas as pd
import numpy as np
import ingest
import networkx as nx



# these classes form a part of the general SV format, TODO move them into dedicated file once the format is finalised
# A position is specified by the organism, chromosome, haplotype and base position
# A range takes a start and an end position. If the end < start, the range is defined as inverted
#TODO use proper cython types, e.g. char for haplo
class Position:
    def __init__(self, org:str, chr:int, haplo:str, pos: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.pos = pos

    def __repr__(self):
        return f"Position({self.org}, {self.chr}, {self.haplo}, {self.pos})"

class Range:
    def __init__(self, org:str, chr:int, haplo:str, start: int, end: int):
        self.org = org
        self.chr = chr
        self.haplo = haplo
        self.start = start
        self.end = end

    def __repr__(self):
        return f"Range({self.org}, {self.chr}, {self.haplo}, {self.start}, {self.end})"

# notes
#   - start with syri output, add alignment info later if needed
#   - identify syntenic regions to ref
#       - find regions sharing similar locations on the A genome, matching on aStart?
#       - sorted join type algorithm?
#       - lifting pansyntenic regions? if lifting, what topology?
#   - output in file/plot in some clever way
#       - maybe as chromosome "bands"?
#       - maybe as % syntenic to each query? => clustering/rooted tree
#       - use plotsr?
#TODO handle positions in non-reference
#TODO output format
#TODO troubleshoot dropped regions
#TO/DO be more lenient?
#TO/DO handle incorrectly mapped chromosomes (use mapping from syri output)

# refactor out the reading in part, TODO move to ingest and interface to a fileformat-agnostic method once file format finalised
def extract_syn_regions(fins, ref="a"):

    # columns to look for as start/end positions
    reforg = "ref"
    refchr = ref + "chr"
    refhaplo = "NaN"
    refstart = ref + "start"
    refend = ref + "end"

    qry = 'b' if ref == 'a' else 'a' # these seem to be the only two values
    qrychr = qry + "chr"
    qryhaplo = "NaN"
    qrystart = qry + "start"
    qryend = qry + "end"

    syns = []
    for fin in fins:
        qryorg = fin.split('/')[-1]
        buf = []
        raw, chr_mapping = ingest.readsyriout(fin)
        raw = raw.loc[raw['type']=='SYN']
        # if implementing filtering later, filter here

        for row in raw.iterrows():
            row = row[1]
            buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
                Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
                ])

        syns.append(pd.DataFrame(data=buf, columns=["ref", qryorg]))

    return syns


def find_pansyn(fins, ref="a"):
    """
    Finds pansyntenic regions by finding the overlap between all syntenic regions in the input files.
    Seems to be very conservative.
    :param: a list of filenames of SyRI output to read in, as well as optionally specifying which sequence is the reference (default 'a').
    :return: a pandas dataframe containing the chromosome, start and end positions of the pansyntenic regions on the reference chromosome, as determined by an overlapping intersection
    """


    syns = extract_syn_regions(fins, ref=ref)
    print(syns)

    # take two dataframes of syntenic regions and compute the flexible intersection
    def intersect_syns(left, right):
        #return pd.merge(left, right) # equivalent to old strict intersection
        
        # they should already be sorted -- required for the rest of the algorithm
        #left.sort_values([refchr, refstart, refend])
        #right.sort_values([refchr, refstart, refend])

        ret = pd.DataFrame(columns = left.columns) # or delete from left? might be faster
        #ret = []

        if len(right) == 0:
            raise ValueError("right is empty!")
        if len(left) == 0:
            raise ValueError("left is empty!")

        riter = right.iterrows()
        rsyn = next(riter)[1]
        liter = left.iterrows()
        lsyn = next(liter)[1]
        while True:
            try: # python iterators suck, so this loop is entirely try-catch'ed
                # this may be very slow, TO/DO benchmark

                # this uses lexicalic comparisons on strings and is most likely fairly slow
                # TO/DO possible improvement: store chromosomes as unsigned byte (cython)
                if rsyn[refchr] > lsyn[refchr]:
                    lsyn = next(liter)[1]
                    continue
                if lsyn[refchr] > rsyn[refchr]:
                    rsyn = next(riter)[1]
                    continue
                
                # determine overlap region
                ovstart = max(rsyn[refstart], lsyn[refstart])
                ovend = min(rsyn[refend], lsyn[refend])

                if ovstart < ovend: # there is valid overlap
                    ret = pd.concat([ret, pd.DataFrame(data=[[lsyn[refchr], ovstart, ovend]], columns=left.columns)])

                # ratchet by dropping the segment with a smaller end
                if lsyn[refend] > rsyn[refend]: # left is after right
                    rsyn = next(riter)[1]
                else: # right is after left
                    lsyn = next(liter)[1]
            except StopIteration: # nothing more to match
                break

        del riter
        del liter
        return ret


    pansyns = functools.reduce(intersect_syns, syns)

    return pansyns


# use networkx for the graph algorithm backend
# how to encode a region? by organism, coordinates & chromosome?
# how to recognize a region? match over organism, find start position/overlap within tolerances?
# => n^2 algorithm probably
# maybe store nodes sorted? networkx uses a dict by default, but can read from lists
# => do it in a (sorted?) list, then use a similar algorithm to the sequence-based code above to identify syntenic regions?
def exact_pansyn(fins):
    """

    :param:
    :return:
    """
    # helper function to extract only large, syntenic regions suitable for analysis from a DF
    def extract_syn_regions(df):
        return df.loc[df['type']=='SYN']

    syns = [extract_syn_regions(ingest.readsyriout(fin)[0]) for fin in fins]


# make it possible to use this file as test
if __name__ == "__main__": # testing
    import sys
    df = find_pansyn(sys.argv[1:])
    print(df)
    print("regions:", len(df))
    print("total lengths:", sum(map(lambda x: x[1][2]-x[1][1],df.iterrows())))
    #df = exact_pansyn(sys.argv[1:])
    #print(df)
    #print("regions:", len(df))
