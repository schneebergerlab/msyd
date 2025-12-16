#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
from collections import defaultdict, deque # or use cpp vector/custom?
from multiprocessing import Pool

import pandas as pd

from msyd.multisyn import Multisyn, Private
from msyd.coords import Range, Panco
from msyd.seq import SeqHandler

import msyd.util as util
import msyd.intersection as intersection

import logging
logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)

# if the distance between two syns is more than this,
# add a private node in between (or do not annotate the link if add_private is not passed)
cdef int MIN_PRIV_THRESH = intersection.get_min_syn_thresh()

cdef class Node:#(Multisyn):
    """
    Internal graph representation for msyd's synteny graph.
    """
    cdef:
        public object index # PanCo
        public str seq
        public object msyn
        public dict[str, Node] post
        public dict[str, Node] prev

    def __init__(self, index, msyn, seqh=None):
        self.index = index
        self.seq = None
        self.msyn = msyn
        self.post = dict()
        self.prev = dict()
        if seqh: # TODO figure out how to provide option for consensus seq etc.
            self.add_sequence(seqh)

    def add_sequence(self, seqh: SeqHandler):
        self.seq = seqh.get_rep_seq(self.msyn)

    def to_gfa1(self, tag_orgs_s=set(), tag_orgs_l=set(), rgfa_tags=True):
        return self.gfa1_s(tag_orgs_s=tag_orgs_s, rgfa_tags=rgfa_tags) + "\n".join(self.gfa1_pre_l(tag_orgs_l=tag_orgs_l))
    
    def gfa1_s(self, tag_orgs_s=set(), rgfa_tags=True):
        ret = f"S\t{self.index}\t*" 
        if tag_orgs_s:
            orgs = tag_orgs_s.union(self.msyn.get_organisms())
            if orgs:
                ret += f"\tSO:Z:{','.join(orgs)}" # originally used ' '
        if rgfa_tags:
            ret += f"\tSN:Z:{self.msyn.ref.org}\tSO:i:{self.msyn.ref.start}\tSR:i:1"
        return ret

    # copy to have post, if necessary later
    def gfa1_pre_l(self, tag_orgs_l=set()):
        # collate all orgs for every previous node
        # allows collapsing all syntenic orgs into one L line
        prevnodes = defaultdict(set)
        for org, node in self.prev.items():
            if not node.is_terminal():
                prevnodes[node].add(org)
        # iterates over all previous nodes, adds the tagged ones as an annotation (if any are tagged)
        return [(f"L\t{node.index}\t+\t{self.index}\t+\t{'*' if not self.seq else self.seq}" if not tag_orgs_l.union(orgs)
                 else f"L\t{node.index}\t+\t{self.index}\t+\t{'*' if not self.seq else self.seq}\tLO:Z:{' '.join(tag_orgs_l.union(orgs))}") for node, orgs in prevnodes.items()]

    #def __hash__(self):
    #    return self.index.__hash__()
    
    def __repr__(self):
        return f"Node({self.index}, {self.msyn}, prev:{list(self.prev.keys())}, post:{list(self.post.keys())})"
    
    def is_terminal(self):
        return self.msyn is None


def make_graphs_chrdict(msyndict, seqh=None, add_private=True, ncores=1):
    """
    Calls make_graph to compute a graph representation from a dictionary containing Multisyn lists indexed by chromosome.
    The graph will be indexed by chromosome and returned in a topological ordering.
    """
    graphs_call = functools.partial(make_graph, seqh=seqh, add_private=add_private)
    ncores = min(len(msyndict), ncores)

    if ncores > 1:
        with Pool(ncores) as pool:
            return dict(pool.map(graphs_call, [msyndict[chrom] for chrom in msyndict]))
    else:
        return dict(map(graphs_call, [msyndict[chrom] for chrom in msyndict]))

cpdef make_graph(msyns, seqh=None, add_private=True):
    """
    Computes a graph representation from a list of Multisyns.
    Coresyns are used to contrain the graph to a single shared node.
    Between two coresyns, nodes are topologically ordered, and their neighbourhood is then reconstructed by tracing along the ordering.
    Returns a DataFrame of Nodes, by default maintaining their topological sorting.
    """
    ret = deque()
    chrom = msyns.iat[0, 0].ref.chr
    # TODO find nice way to figure out if incrementing core or mera counter; or do panco imputation after graph construction
    #panco = Panco(chrom, 0, 0)
    curdict = dict()#defaultdict(lambda: None)

    logger.info(f"Starting graph construction on {chrom}, containing {msyns.size} Msyns")
    msyns.sort_values(by=[0], inplace=True) # note: if too inefficient, fetch regions between coresyns and call subfunction to sort, or implement insertion sort
    logger.info(f"Finished top. sorting on {chrom}")
    util.validate_top_sort(msyns)

    # store start and end, add as virtual nodes later
    starting = Node(None)
    ending = Node(None)

    logger.info(f"Starting graph construction on {chrom}")
    for _, msyn in msyns.iterrows():
        msyn = msyn[0]
        node = Node(msyn, seqh=seqh)
        # add links to predecessors per organism
        for org, rng in msyn.iter_orgs_ranges(): #[(msyn.ref.org, msyn.ref)] + list(msyn.ranges_dict.items()): # how to handle ref?
            if org in curdict: # default case
                curprev = curdict[org]
                curprevrng = curprev.msyn.ranges_dict[org] if org in curprev.msyn.ranges_dict else curprev.msyn.ref # has to be on ref if it isn't in ranges_dict
                # make sure this assumption is valid
                assert org == curprevrng.org

                # check distance, add direct link or private region
                if rng.start - curprevrng.end < MIN_PRIV_THRESH:
                    # add link and backlink
                    node.prev[org] = curprev
                    curprev.post[org] = node
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
                #node.prev['_start'] = starting # to make the graph traversable
        # done with looping over orgs
        ret.append(node)

    # log current state as ending nodes
    ending.pre = curdict
    #for node in curdict.values():
    #    node.post['_end'] = ending

    ret.appendleft(ending) # second pos
    ret.appendleft(starting) # first pos

    return (chrom, pd.DataFrame(data=list(ret)))

cpdef trace_org(begin, org, forward=True):
    """
    Traces an organisms path through the graph.
    Starts at the node begin, which has to contain org.
    Traces in either forward (default) or backward direction (taking the post/prev dict each time) depending on the parameter passed to forward.
    Returns a List containing the nodes org traverses through in order.
    """
    ret = deque()
    cur = begin

    while True: # graph is a DAG, no need to worry about cycles
        if cur.msyn: # to not append start/end node
            # make sure we don't mistraverse
            #assert org in cur.msyn.get_organisms()
            ret.append(cur)

        # continue traversal
        iterdict = cur.post if forward else cur.prev
        if iterdict:
            cur = iterdict[org]
        else:
            logger.info(f"Finished traversing on {org} at {cur}")
            break

    return list(ret)

cpdef trace_bidirectional(begin, org):
    # should also work if called on start/end, indexing defaults to []
    return trace_org(begin, org, forward=False) + trace_org(begin, org, forward=True)[1:]
