# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import pandas as pd
import pysam

import io
import os
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

    @classmethod
    def from_fasta_tsv(cls, fin):
        """
        Instantiates a `SeqHandler` from a TSV containing the organism name in the first column and their correpsonding FASTA file in the second column.
        """
        # check input is file-like/can be opened
        if isinstance(fin, (str, os.PathLike)):
            fin = open(fin, 'rt')
        elif not isinstance(fin, io.TextIOBase):
            raise ValueError(f"{fin} is not a path-like or file-like object!")

        _backing = dict()

        # parse file line by line
        for line in fin:
            if line[0] == '#' or line.strip() == '':
                continue
            cells = line.strip().split('#')[0].split('\t')
            if len(cells) != 2:
                logger.error(f"invalid entry in {fin.name}: '{line}' does not contain two columns. Skipping!")
                continue

            org = cells[0].strip()
            fasta = cells[1].strip()
            if not os.path.isfile(fasta):
                raise FileNotFoundError(f"Cannot find file at {fasta}. Double-check the input TSV. Exiting.")
            _backing[org] = pysam.FastaFile(fasta)

        return cls(_backing)

    def get_len(self, org, chrom):
        return self._backing[org].get_reference_length(chrom)

    def get_len_dict(self, chrom):
        return {org: self._backing[org].get_reference_length(chrom) for org in self._backing}

    def get_fn_dict(self):
        """
        Returns this SeqHandler in fasta dict format.
        Can be passed to from_fasta_dict to rehydrate the SeqHandler.
        """
        return {org: self._backing[org].filename for org in self._backing}

    # support pickling, for use with multiprocessing
    def __getstate__(self):
        return self.get_fn_dict()

    def __setstate__(self, state):
        # call constructor and steal its _backing, to ensure consistent call
        self._backing = self.from_fasta_dict(state)._backing 
        #{org: pysam.FastaFile(fasta) for org, fasta in state.items()}

    def get_range(self, rng: Range, margin=0):
        """
        Retrieves the sequence corresponding to the positions within this range on the organism it is described on.
        """
        assert rng.org in self._backing
        fa = self._backing[rng.org]
        end = min(fa.get_reference_length(rng.chr), rng.end + margin + 1)
        
        #logger.info(f"Fetching {rng} with {margin} bp margins. End at {end}.")

        return fa.fetch(region=rng.chr, start=max(0, rng.start-margin), end=end).upper()

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

