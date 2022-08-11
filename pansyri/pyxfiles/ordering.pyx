#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import math

import pansyri.ingest as ingest
import pansyri.util as util

def order(syns, alns):
    """Convenience function performing a full ordering imputation given syns/alns extracted from a.tsv
    Mostly meant to be called from main.py.
    :param syns, alns: list of syn/aln files, as output by parse_tsv in util
    :returns: the calculated ordering of the organisms as a `list` of strings
    :rtype: List[str]
    """
    df = util.crosssyn_from_lists(syns, alns, SYNAL=False, cores=6)
    print("INFO: got crossyn df")
    orgs = util.get_orgs_from_df(df)
    print("INFO: got orgs from crossyn df")
    return order_greedy(orgs, df)

def syn_score(cur, org, df):
    """Defines a similarity score from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    # using length on org, but this shouldn't matter too much
    return sum(map(lambda x: len(x.ranges_dict[org]),filter(lambda x: cur in x.ranges_dict and org in x.ranges_dict, map(lambda x: x[1][0], df.iterrows()))))
    # Ideas for future better synscores:
    # – filter for degree (maybe lower and upper bound? likely the medium degree crosssyntenic regions most informative
    # – have some non-linearity in the length processing, e.g. square it
    # – maybe combine somehow with variant-based scoring?

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

    cur = list(orgs)[0] # arbitrarily choose first organism
    order = [cur]
    orgs.remove(cur)

    while orgs:
        print("INFO: currently selected:", orgs)
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


#def len_correct(score_fn):
#    """Higher-order function returning a length-corrected scoring function given an uncorrected score.
#    """
#    def corrected(syri):
#        gen1_len = sum(map(lambda x: len(x), ingest.readfasta(gen1).values()))
#        gen2_len = sum(map(lambda x: len(x), ingest.readfasta(gen2).values()))
#        l_eff = (gen1_len + gen2_len)/2
#        return score_fn(syri)/l_eff
#    return corrected
