#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3
import math

import logging

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import complete, dendrogram, leaves_list

import pansyri.util as util

logger = logging.getLogger(__name__)

def order(syns, alns, chr=None):
    """Convenience function performing a full ordering imputation given syns/alns extracted from a.tsv
    Mostly meant to be called from main.py.
    It also demonstrates the way to compute an ordering using the new procedure:
    First, the crosssyntenic regions are extracted from the SyRI and alignment files, supplied e.g. by reading from a tsv file.
    Then, the organism names are extracted from the Crosssynteny DF (this might be moved to the ingest module in the future).
    Finally, the optimizing algorithm is called with the appropriate scoring function, in this case the greedy algorithm with the default syn_score.
    :param syns, alns: list of syn/aln files, as output by parse_tsv in util
    :param chr: The chromosome to restrict the ordering to. If `None` (default), the ordering is performed across the entire genome. The filtering is performed according to chromosome on the reference.
    :type chr: `str`
    :returns: the calculated ordering of the organisms as a `list` of strings
    :rtype: List[str]
    """
    # disabled using SYNAL for now, runtime fairly slow on the dell nodes, not quite sure why
    # might even be the better idea, we want to capture large-scale synteny anyway
    #df = util.crosssyn_from_lists(syns, alns, SYNAL=False, cores=6)

    # hack together union of coresyns, not elegant or fast but should get the job done. Runtime in Θ(n^2)!!
    import itertools
    df = pd.concat([util.coresyn_from_lists(list(synpair), list(alnpair), SYNAL=False) for synpair, alnpair in zip(itertools.product(syns, syns), itertools.product(alns, alns)) if synpair[0] != synpair[1]])

    #logger.debug(df.head(100).to_string())
    logger.info("got crossyn df")

    if chr is not None:
        df = util.filter_multisyn_df_chr(df, chr)
        # if filtering to a range of interest, call filter_multisyn_df instead like this:
        # df = util.filter_multisyn_df(df, Range(None, <chr>, 'NaN', <start>, <end>))
        logger.info("Filtered Crossyn DF")

    # optionally adjust the organism list
    # orgs = util.get_orgs_from_df(df)[:4]
    # make the call to the optimizer, using the default/only scoring function right now:
    return order_complete(df, orgs=None, score_fn=syn_score)

def syn_score(cur, org, df):
    """Defines a similarity score from syri output.
    Currently just uses the sum of the lengths of all syntenic regions without correcting for genome size.
    """
    if org == 'ref':    
        return sum(map(lambda x: len(x.ranges_dict[cur]), filter(lambda x: cur in x.ranges_dict, map(lambda x: x[1][0], df.iterrows()))))

    if cur == 'ref':    
        return sum(map(lambda x: len(x.ranges_dict[org]), filter(lambda x: org in x.ranges_dict, map(lambda x: x[1][0], df.iterrows()))))

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

def order_complete(df, orgs=None, score_fn=syn_score, maximize=True, ref=True):
    if orgs is None:
        logger.info("getting orgs from crossyn df as none were supplied")
        orgs = util.get_orgs_from_df(df)
    if ref is True:
        orgs.add('ref') # include reference
    orgs = sorted(set(orgs))
    n = len(orgs)

    # construct and fill the distance matrix
    distmat = np.zeros([n, n])
    for x in range(n):
        logger.debug(f"Distance matrix calculation {x}/{n}")
        for y in range(x+1, n):            
            distmat[x][y] = score_fn(orgs[x], orgs[y], df)

    # make the scipy call
    Z = complete(distmat)
    order = [orgs[i] for i in leaves_list(Z)]
    return order
    


def order_greedy(df, orgs=None, score_fn=syn_score, maximize=True, ref=True):
    """A simple, greedy algorithm ordering a list of organisms while trying to maximize (or minimize, depending on `maximize`) the similarity score between each organism and the next one.
    The first organism is chosen at random.

    :param df: a `DataFrame` of `Pansyn` objects, as produced by `find_multisyn(detect_crosssyn=True)`.
    :param orgs: a sequence of organism names. If `None`, all the org names in `df` are used.
    :type orgs: Set[str]

    :param score_fn: similarity score to use, by default `syn_score`. The scoring function needs to accept the filename of a syri.out file as output.
    :param_type score_fn: a function mapping a file path to a numerical score
    :param maximize: Boolean toggling whether to minimize/maximize the score (similarity vs dissimilarity scoring). Defaults to True.
    :param ref: `bool` controlling whether to include the reference in the ordering
    :param_type ref: `bool`
    :returns: a list containing the elements of orgs ordered according to the algorithm.
    :rtype: List[str]
    """
    if orgs is None:
        logger.info("getting orgs from crossyn df as none were supplied")
        orgs = util.get_orgs_from_df(df)
    if ref is True:
        orgs.add('ref') # include reference
    orgs = set(orgs)

    cur = list(orgs)[0] # arbitrarily choose first organism
    print(orgs)
    order = [cur]
    orgs.remove(cur)

    while orgs:
        # find the next organism with maximal similarity score to this one
        ext_score = -math.inf if maximize else math.inf
        for org in orgs:
            score = score_fn(cur, org, df)
            logger.debug(f"{cur}, {org}: {score}")
            if maximize and score > ext_score:
                ext_score = score
                cur = org
            elif not maximize and score < ext_score:
                ext_score = score
                cur = org

        logger.debug(f"selected: {order}, score: {ext_score}")
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
