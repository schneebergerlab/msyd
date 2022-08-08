#!/usr/bin/python3
# maybe cythonize later, probably not worth it though
import pansyri.util as util
#import math

def syn_score(syri):
    """
    Defines a distance metric from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    syns = util.extract_regions(syri, ann='SYN')
    return sum(map(lambda x: len(x[1][0]), syns.iterrows()))

def order_plotsr_greedy(orgs, score_fn=syn_score, filename_mapper=lambda x, y: x+'_'+y+"syri.out"):
    """
    A simple, greedy algorithm ordering a list of organisms while trying to maximize the similarity score between each organism and the next one.
    :params:
        orgs is a sequence of organism/filenames
        score sets the similarity score to use, by default `syn_score`. The scoring function needs to accept the filename of a syri.out file as output.
        filename_mapper turns the names of two organisms into the filename of the corresponding syri output.
    :returns: a list containing the elements of orgs ordered according to the algorithm.
    """
    orgs = set(orgs)

    cur = list(orgs)[0] # arbitrarily choose start point
    order = [cur]
    orgs.remove(cur)

    while orgs:
        # find the next organism with maximal similarity score to this one
        max_score = 0# math.inf
        for org in orgs:
            score = score_fn(filename_mapper(cur, org))
            if score > max_score:
                max_score = score
                cur = org

        orgs.remove(cur)
        order.append(cur)

    return order

