#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

from msyd.multisyn import Multisyn, Private
from msyd.coords import Range, Panco
from msyd.utils import *

import functools
from collections import defaultdict, deque # or use cpp vector/custom?
from multiprocessing import Pool

import pandas as pd

import logging
logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)

# if the distance between two syns is more than this,
# add a private node in between (or do not annotate the link if add_private is not passed)
cdef int MIN_PRIV_THRESH = intersection.get_min_syn_thresh()

cdef class Node:
    """
    Internal graph representation for msyd's synteny graph.
    """
    #cdef:
    #    public Multisyn msyn
    #    public dict[str, Node] post
    #    public dict[str, Node] prev
    #    public Panco index

    def __init__(self, msyn):
        self.msyn = msyn
        self.post = dict()
        self.prev = dict()

    cdef to_gfa1(self):
        #TODO implement serialization as one S line and L lines to successors
        #TODO think about adding path lines for every org at the end in another function
        return self.gfa1_s() +"\n" + "\n".join(self.gfa1_pre_l())
    
    cdef gfa1_s(self, tag_orgs_s=set(), rgfa_tags=True):
        #TODO fetch sequence somehow
        ret = f"S\t{self.index}\t*" 
        if tag_orgs_s:
            orgs = tag_orgs_s.union(self.msyn.get_organisms())
            if orgs:
                ret += f"\tSO:Z:{' '.join(orgs)}"
        if rgfa_tags:
            ret += f"SN:Z:{self.msyn.ref.org}\tSO:i:{self.msyn.ref.start}\tSR:i:1"
        return ret

    # copy to have post, if necessary later
    cdef gfa1_pre_l(self, tag_orgs_s=set()):
        # collate all orgs for every previous node
        # allows collapsing all syntenic orgs into one L line
        prevnodes = defaultdict(set)
        for org, node in prev.items():
            if node in prevnodes:
                prevnodes[node].append(org)
        # iterates over all previous nodes, adds the tagged ones as an annotation (if any are tagged)
        return [(f"L\t{node.index}\t+\t{self.index}\t+\t*" if not tag_orgs_s.union(orgs)
                 else f"L\t{node.index}\t+\t{self.index}\t+\t*\tLO:Z:{' '.join(tag_orgs_l.union(orgs))}") for node, orgs in prev.items()]

    def __hash__(self):
        return self.index


def make_graphs_chrdict(msyndict, add_private=True, cores=1):
    """
    Calls make_graph to compute a graph representation from a dictionary containing Multisyn lists indexed by chromosome.
    The graph will be indexed by chromosome and returned in a topological ordering.
    """
    graphs_call = functools.partial(make_graph, add_private=add_private)
    cores = min(len(msyndict), cores)

    if cores > 1:
        with Pool(cores) as pool:
            return dict(pool.map(graphs_call, [msyndict[chrom] for chrom in syndict]))
    else:
        return dict(map(graphs_call, [msyndict[chrom] for chrom in syndict]))

cpdef make_graph(msyns, add_private=True):
    """
    Computes a graph representation from a list of Multisyns.
    Coresyns are used to contrain the graph to a single shared node.
    Between two coresyns, nodes are topologically ordered, and their neighbourhood is then reconstructed by tracing along the ordering.
    Returns a DataFrame of Nodes, by default maintaining their topological sorting.
    """
    logger.info(f"Starting graph construction on {msyns.iloc[0].ref.chr}, containing {msyns.size} Msyns")
    ret = deque()
    chrom = msyns.iloc[0].ref.chrom
    # TODO find nice way to figure out if incrementing core or mera counter; or do panco imputation after graph construction
    #panco = Panco(chrom, 0, 0)
    curdict = dict()#defaultdict(lambda: None)
    msyns.sort_values(inplace=True) # note: if too inefficient, fetch regions between coresyns and call subfunction to sort, or implement insertion sort
    logger.info(f"Finished top. sorting on {chrom}")
    util.validate_top_sort(msyns)

    # store these for now, not sure how best to handle
    # maybe annotate as virtual node w/o msyn? => easier for path tracing
    starting = Node(None)
    ending = Node(None)

    logger.info(f"Starting graph construction on {chrom}")
    for _, msyn in msyns.iterrows():
        msyn = msyn[0]
        node = Node(msyn)
        # add links to predecessors per organism
        for org, rng in [(msyn.ref.org, msyn.ref)] + msyn.ranges_dict.items(): # how to handle ref?
            if org in curdict: # default case
                curprev = curdict[org]
                curprevrng = curprev.ranges_dict[org] if org in curprev.ranges_dict else curprev.ref # has to be on ref if it isn't in ranges_dict

                # check distance, add direct link or private region
                if rng.start - curprevrng.end < MIN_PRIV_THRESH:
                    # add link and backlink
                    node.prev[org] = curdict[org]
                    curprev.post = node
                else: # add private region
                    privnode = Node(Private(Range(org, chrom, curprevrng.end + 1, rng.start -1)))

                    # add two back/frontlinks
                    curprev.post[org] = privnode
                    privnode.prev[org] = curprev
                    privnode.post[org] = node
                    node.prev[org] = privnode

                    # add private node to node list
                    ret.append(privnode)
                # change curdict to this node
                curdict[org] = node
            else: # init case
                starting.post[org] = node
                curdict[org] = node
        # done with looping over orgs
        ret.append(node)

    # log current state as ending nodes
    ending.pre = curdict

    return pd.DataFrame(data=[starting, ending] + ret)

cpdef trace_org(begin, org, forward=True):
    """
    Traces an organisms path through the graph.
    Starts at the node begin, which has to contain org.
    Traces in either forward (default) or backward direction (taking the post/prev dict each time) depending on the parameter passed to forward.
    Returns a List containing the ordered nodes org traverses through.
    """
    if not org in begin.msyn.get_organisms():
        raise ValueError(f"{org} not found in starting Node {begin}!")
    return []
