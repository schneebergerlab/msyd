#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import math

#TODO do all qry to reference as rest of pansyri
# use output from find_multisyn
# can do either pairwise for each or extract from whole output
# pairwise might be faster? but maybe irrelevant given parallelization
# full has advantage of selecting for degrees, and later multigenomic SVs

# to refactor:
# – use same input as main, parse the tsv using util.parse_tsv
# – remove size correction for now => no more need for genomes
    # => probably not all that relevant if using subset from range anyway, should have equal lengths
    # possibly readd with length from ranges later?
# – filter for range, matching on chromosome
# for now just do that here using DataFrame.loc[]


import pansyri.ingest as ingest
import pansyri.util as util

def order(syns, alns):
    df = util.crosssyn_from_lists(syns, alns, cores=6)
    orgs = util.get_orgs_from_df(df)
    return order_greedy(orgs, df)

def syn_score(cur, org, df):
    """Defines a similarity score from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    # using length on org, but this shouldn't matter too much
    return sum(map(lambda x: len(x.ranges_dict[org]),filter(lambda x: cur in x.ranges_dict and org in x.ranges_dict, map(lambda x: x[1][0], df.iterrows()))))

#def len_correct(score_fn):
#    """Higher-order function returning a length-corrected scoring function given an uncorrected score.
#    """
#    def corrected(syri):
#        gen1_len = sum(map(lambda x: len(x), ingest.readfasta(gen1).values()))
#        gen2_len = sum(map(lambda x: len(x), ingest.readfasta(gen2).values()))
#        l_eff = (gen1_len + gen2_len)/2
#        return score_fn(syri)/l_eff
#    return corrected

# temporarily commented out, later reincorporate from proper multigenomic varcalling
# needs Equality between SVs solved, redone syri imports
#def sv_score(cur, org, df):
#    """Defines a dissimilarity score by summing up the length of all regions annotated to be a structural variant, excluding duplications.
#    """
#    return sum(map(lambda x: len(x[1][0]), ingest.extract_syri_regions(syri, anns=['INV', 'TRANS', 'DEL', 'INS']).iterrows()))


def order_greedy(orgs, df, score_fn=syn_score, maximize=True):
    """A simple, greedy algorithm ordering a list of organisms while trying to maximize (or minimize, depending on `maximize`) the similarity score between each organism and the next one.
    The first organism is chosen at random.

    :param orgs: a sequence of organism names.
    :param df: a DataFrame of `Pansyn` objects, as produced by `find_multisyn(detect_crosssyn=True)`.
    :param score_fn: similarity score to use, by default `syn_score`. The scoring function needs to accept the filename of a syri.out file as output.
    :param_type score_fn: a function mapping a file path to a numerical score
    :param maximize: Boolean toggling whether to minimize/maximize the score (similarity vs dissimilarity scoring). Defaults to True.
    :returns: a list containing the elements of orgs ordered according to the algorithm.
    :rtype: List[str]
    """
    orgs = set(orgs)

    cur = list(orgs)[0] # arbitrarily choose first organism
    order = [cur]
    orgs.remove(cur)

    while orgs:
        # find the next organism with maximal similarity score to this one
        ext_score = -math.inf if maximize else math.inf
        for org in orgs:
            score = score_fn(cur, org, df)
            if maximize and score > ext_score:
                ext_score = score
                cur = org
            elif not maximize and score < ext_score:
                ext_score = score
                cur = org

        orgs.remove(cur)
        order.append(cur)

    return order

