#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

class SNV:
    """
    """
    def __init__(self, ref, qry, refseq=None, qryseq=None):
        self.ref = ref
        self.qry = qry
        self.refseq = refseq
        self.qryseq = qryseq
