#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import sys
from collections import defaultdict, deque # or use cpp vector/custom?

from typing import TypeAlias, List, Dict

import msyd.util as util

import logging
logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)


Org: TypeAlias = str

class OrgContainer:
    #public List[Org] orgs
    #public int n
    #public Dict[Org, Int] posdict

    def __init__(self):
        self.orgs = []
        self.n = 0 # keep track of how big orgs is, to save having to call len() often
        self.posdict = {} # track position of each ref; necessary?

    def get_n(self) -> int:
        """
        Returns the no. of organisms in this population.
        """
        return self.n

    def get_orgs(self) -> List[Org]:
        """
        Returns the names of the organisms in this population, in order.
        Meant to allow iteration across all organisms, and may change to a generator at any point.
        To work with specific organism from this container, use `get_pos` or `get_at`.
        """
        return self.orgs

    def get_at(self, pos: int) -> Org:
        return self.orgs[pos]

    def get_ind(self, org: Org) -> int:
        """
        Get the index of a given org in this OrgContainer.
        """
        if not org in self.orgs:
            logger.error(f"OrgContainer queried with missing org: {org}")
            raise ValueError(f"OrgContainer queried with missing org: {org}")
        return self.posdict[org]

    def add_org(self, org:str):
        org = sys.intern(org) # make sure this goes to the string literal pool
        self.posdict[org] = self.n
        self.orgs.append(org)
        self.n += 1

    def sort(self):
        """
        Sorts this OrgContainer alphabetically.
        Mutates self!
        """
        self.orgs = sorted(self.orgs)
        self.posdict = {key:pos for pos, key in enumerate(self.orgs)}

    # its convenient, but may confuse Cython code, so don't implement
    #def __getitem__(self, key):
    #    if isinstance(key, int):
    #        return self.get_pos(key)
    #    else:
    #        return self.posdict[key]

    def __contains__(self, org:Org):
        return self.posdict.__contains__(org) # should be faster to ask the hashset than the list

    @classmethod
    def from_list(cls, lst):
        ret = cls()
        for org in lst:
            ret.add_org(org)
        return ret






