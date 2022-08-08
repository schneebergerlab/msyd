#!/usr/bin/python3
# maybe cythonize later, probably not worth it though
import pansyri.util as util
import math

def syn_score(syri):
    """
    Defines a distance metric from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    syns = util.extract_regions(syri, ann='SYN')
    return sum(map(lambda x: len(x[1][0]), syns.iterrows()))

def 

def order_greedy(orgs, score_fn=syn_score, filename_mapper=lambda x, y: x+'_'+y+"syri.out", maximize=True):
    """A simple, greedy algorithm ordering a list of organisms while trying to maximize (or minimize, depending on `maximize`) the similarity score between each organism and the next one.
    :param orgs: a sequence of organism/filenames
    :param score_fn: similarity score to use, by default `syn_score`. The scoring function needs to accept the filename of a syri.out file as output.
    :param_type score_fn: a function mapping a file path to a numerical score
    :param filename_mapper: a function mapping two genome ids to a syri file comparing the two. By default set up to work with the `../data/ampril` dataset.
    :param_type filename_mapper: lambda str, str: str
    :param maximize: Boolean toggling whether to minimize/maximize the score (similarity vs dissimilarity scoring). Defaults to True.
    :returns: a list containing the elements of orgs ordered according to the algorithm.
    :rtype: List[str]
    """
    orgs = set(orgs)

    cur = list(orgs)[0] # arbitrarily choose start point
    order = [cur]
    orgs.remove(cur)

    while orgs:
        # find the next organism with maximal similarity score to this one
        ext_score = 0 if maximize else math.inf
        for org in orgs:
            score = score_fn(filename_mapper(cur, org))
            if maximize and score > ext_score:
                ext_score = score
                cur = org
            elif not maximize and score < ext_score:
                ext_score = score
                cur = org

        orgs.remove(cur)
        order.append(cur)

    return order

