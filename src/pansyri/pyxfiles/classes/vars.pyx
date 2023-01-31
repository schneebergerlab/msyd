#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pansyri.util as util

logger = util.CustomFormatter.getlogger(__name__)

class SNV:
    """A class storing a SNV.
    TODO how to expand to multigenomics? List[SNV]? SNV{qry=List/dict[Position]}? how to store multiple different alleles?
    :param ref:
    :type ref: `Position`
    :param qry:
    :type qry: `Position`
    :param refseq:
    :type refseq: `str`
    :param qryseq:
    :type qryseq: `str`
    """
    def __init__(self, ref, qry, refseq=None, qryseq=None):
        self.ref = ref
        self.qry = qry
        self.refseq = refseq
        self.qryseq = qryseq


class Inv:
    pass

class Dup:
    pass

class Ins:
    pass

class Del:
    pass
