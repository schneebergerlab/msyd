# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import pysam

from typing import Dict

from msyd.coords import Range
from msyd.multisyn import Multisyn
import msyd.util as util

import logging
logger = util.CustomFormatter.getlogger(__name__)
logger.setLevel(logging.INFO)

class SeqHandler:
    """
    Manages genomic sequences of multiple organisms.
    Uses pysam's FastaFile API, which in theory should mean that not everything needs to be kept in memory at once.
    """
    # Notes
    # should genomic sequences be treated as byte arrays or C[++] strings instead? => avoid string literal pool clogging

    def __init__(self, _backing: Dict[str, pysam.FastaFile]):
        self._backing = _backing

    @classmethod
    def from_fasta_dict(cls, fasta_dict: Dict[str, str]):
        """
        Instantiates a `SeqHandler` from a dictionary mapping organism names to a FASTA file of their genome.
        """
        _backing = {org: pysam.FastaFile(fasta) for org, fasta in fasta_dict.items()}
        return cls(_backing)

    def get_range(self, rng: Range, margin=0):
        """
        Retrieves the sequence corresponding to the positions within this range on the organism it is described on.
        """
        assert rng.org in self._backing
        fa = self._backing[rng.org]
        end = min(fa.get_reference_length(rng.chr), rng.end + margin + 1)
        return self._backing.fetch(reference=rng.chr, start=max(0, rng.start-margin), end=end).upper()

    def get_all_seq(self, msyn: Multisyn):
        """
        Return a dictionary containing the sequence of `msyn` on all organisms it is found on.
        """
        return {org: self.get_range(rng) for org, rng in msyn.iter_orgs_ranges()}

    def get_rep_seq(self, msyn: Multisyn):
        """
        Get a representative sequence associated for `msyn`.
        Defaults to the sequence on the genome chosen as the common alignment target.
        """
        return self.get_range(msyn.ref)

    def get_consensus_seq(self, msyn: Multisyn):
        """
        Get a consensus sequence for an `msyn`.
        Currently not implemented.
        """
        raise NotImplemented("TODO")

