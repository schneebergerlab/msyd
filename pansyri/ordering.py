#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import math

import pansyri.ingest as ingest

def syn_score(syri, gen1, gen2):
    """Defines a similarity score from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    return sum(map(lambda x: len(x[1][0]), ingest.extract_syri_regions(syri, anns=['SYN']).iterrows()))

def len_correct(score_fn):
    """Higher-order function returning a length-corrected scoring function given an uncorrected score.
    """
    def corrected(syri, gen1, gen2):
        gen1_len = sum(map(lambda x: len(x), ingest.readfasta(gen1).values()))
        gen2_len = sum(map(lambda x: len(x), ingest.readfasta(gen2).values()))
        l_eff = (gen1_len + gen2_len)/2
        return score_fn(syri, gen1, gen2)/l_eff
    return corrected

def sv_score(syri, gen1, gen2):
    """Defines a dissimilarity score by summing up the length of all regions annotated to be a structural variant, excluding duplications.
    """
    return sum(map(lambda x: len(x[1][0]), ingest.extract_syri_regions(syri, anns=['INV', 'TRANS', 'DEL', 'INS']).iterrows()))


def order_greedy(orgs, score_fn=syn_score, gen_mapper=lambda x: x + '.filtered.fa', syri_mapper=lambda x, y: x+'_'+y+"syri.out", maximize=True):
    """A simple, greedy algorithm ordering a list of organisms while trying to maximize (or minimize, depending on `maximize`) the similarity score between each organism and the next one.
    The first organism is chosen at random.
    :param orgs: a sequence of organism/filenames
    :param score_fn: similarity score to use, by default `syn_score`. The scoring function needs to accept the filename of a syri.out file as output.
    :param_type score_fn: a function mapping a file path to a numerical score
    :param gen_mapper: a function mapping a genome id to a fasta file containing the genome corresponding to that ID. By default set up to work with the `../data/ampril` dataset.
    :param_type filename_mapper: lambda str: str
    :param syri_mapper: a function mapping two genome ids to a syri file comparing the two. By default set up to work with the `../data/ampril` dataset.
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
            score = score_fn(syri_mapper(cur, org), gen_mapper(cur), gen_mapper(org))
            if maximize and score > ext_score:
                ext_score = score
                cur = org
            elif not maximize and score < ext_score:
                ext_score = score
                cur = org

        orgs.remove(cur)
        order.append(cur)

    return order

