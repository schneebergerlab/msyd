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
        return self.n

    def get_orgs(self) -> List[Org]:
        return self.orgs

    def get_pos(self, org: Org) -> int:
        if not org in self.orgs:
            logger.error(f"OrgContainer queried with missing org: {org}")
            raise ValueError(f"OrgContainer queried with missing org: {org}")
        return self.posdict[org]

    def add_org(self, org:str):
        org = sys.intern(org) # make sure this goes to the string literal pool
        self.posdict[org] = self.n
        self.orgs.append(org)
        self.n += 1

    @classmethod
    def from_list(cls, lst):
        ret = cls()
        for org in lst:
            ret.add_org(org)
        return ret






