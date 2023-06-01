#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import functools
import multiprocessing
import pandas as pd
from collections import deque
import os
import io
import logging
import re
import tempfile
import random

cdef TMPDIR = None

class CustomFormatter(logging.Formatter):
    '''
    Copied from https://github.com/mnshgl0110/hometools/blob/master/hometools/classes.py
    who copied from
    https://betterstack.com/community/questions/how-to-color-python-logging-output/
    '''

    grey = "\x1b[0;49;90m"
    green = "\x1b[0;49;32m"
    yellow = "\x1b[0;49;93m"
    red = "\x1b[0;49;31m"
    bold_red = "\x1b[0;49;31;21m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s (%(module)s:%(lineno)d)"
    #TODO why is module always main?? maybe try instantiating not using root logger??

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
    

    def getlogger(name):
        logger = logging.getLogger(name.split('.')[-1])
        handler = logging.StreamHandler()
        handler.setFormatter(CustomFormatter())
        logger.addHandler(handler)
        logging.basicConfig(level=logging.INFO)
        logger.propagate = False
        return logger

logger = CustomFormatter.getlogger(__name__)

# copied from https://stackoverflow.com/questions/50878960/parallelize-pythons-reduce-command
# doesn't seem to be very fast?
def parallel_reduce(reduceFunc, l, numCPUs):
    if numCPUs == 1 or len(l) <= 3:
            returnVal = functools.reduce(reduceFunc, l[1:], l[0])
            return returnVal

    parent1, child1 = multiprocessing.Pipe()
    parent2, child2 = multiprocessing.Pipe()
    p1 = multiprocessing.Process(target=parallel_reduce, args=(reduceFunc, l[:len(l) // 2], numCPUs // 2, child1, ) )
    p2 = multiprocessing.Process(target=parallel_reduce, args=(reduceFunc, l[len(l) // 2:], numCPUs // 2 + numCPUs%2, child2, ) )
    p1.start()
    p2.start()
    leftReturn, rightReturn = parent1.recv(), parent2.recv()
    p1.join()
    p2.join()
    returnVal = reduceFunc(leftReturn, rightReturn)
    return returnVal


cdef gettmpfile():
    if not TMPDIR:
        return tempfile.NamedTemporarzFile().name
    else:
        randstr = ''.join(random.choices('abcdefghijklmnopqrstvwxyz', k=6))
        path = f"{TMPDIR}/tmp{randstr}"
        if os.path.isfile(path):
            logger.error(f'Temp file path already exists: {path}'
        open(path, 'w').close() # create the file as empty before returning the path
        return path


#############################################
# everything above this doesn't depend on pasy
#############################################

import pasy.pansyn as pansyn
from pasy.classes.cigar import Cigar
from pasy.classes.coords import Pansyn, Range
import pasy.io



def parse_input_tsv_path(path):
    """DEPRECATED
    Convenience wrapper to call parse_input_tsv with a path.
    """
    with open(path, 'r') as fin:
        return parse_input_tsv(fin)

def parse_input_tsv(fin):
    """
    Takes a tsv file containing the input names/alignments/syri/vcf files and processes it for find_multisyn.
    Anything after a # is ignored. Lines starting with # are skipped.
    By convention, tsv files usually have a header line starting with # (though this is not required).
    :params: `os.PathLike`, `str` or a TextIO object containing the paths of the input alignment and syri files in tsv format.
    Will be consumed by this function!
    :returns: a tuple of three lists containing the organism names, paths of the syri and alignment files, respectively.
    """
    if isinstance(fin, (str, os.PathLike)):
        fin = open(fin, 'rt')
    elif not isinstance(fin, io.TextIOBase):
        raise ValueError(f"{fin} is not a path-like or file-like object!")

    cdef:
        qrynames = deque()
        syris = deque()     # Lists are too slow appending, using deque instead
        alns = deque()
        vcfs = deque()
        fastas = deque()
        str qry = ''
        str syri = ''
        str aln = ''
        str vcf = ''
        str fasta = ''

    for line in fin:
        if line[0] == '#' or line.strip() == '':
            continue

        cells = line.strip().split('#')[0].split('\t')
        if len(cells) < 3:
            logger.error(f"invalid entry in {fin.name}: '{line}' too short. Skipping!")
            continue
        elif len(cells) < 4:
            logger.warning(f"No vcf specified in {fin.name}, trying to find syri vcf. Offending line: {line}")
            qry = cells[0].strip()
            aln = cells[1].strip()
            syri = cells[2].strip()
            vcf = ''.join(syri.split('.')[:-1]) + '.vcf'
        else:
            if len(cells) > 4:
                logger.warning(f"More than five columns in {fin.name}, ignoring anything after fourth column")
            qry = cells[0].strip()
            aln = cells[1].strip()
            syri = cells[2].strip()
            vcf = cells[3].strip()
            fasta = cells[4].strip()

        # Check that the files are accessible
        if not os.path.isfile(aln):
            raise FileNotFoundError(f"Cannot find file at {aln}. Exiting")
        if not os.path.isfile(syri):
            raise FileNotFoundError(f"Cannot find file at {syri}. Exiting")
        if not os.path.isfile(vcf):
            raise FileNotFoundError(f"Cannot find file at {vcf}. Exiting")
        if not os.path.isfile(fasta):
            raise FileNotFoundError(f"Cannot find file at {fasta}. Exiting")

        qrynames.append(qry)
        alns.append(aln)
        syris.append(syri)
        vcfs.append(vcf)
        fastas.append(fasta)

    if len(set(qrynames)) != len(qrynames):
        logger.error(f"Non-unique names in {fin.name}. This will most likely cause issues, proceed with caution!")

    fin.close()
    return (qrynames, syris, alns, vcfs, fasta)


# set of utility funcitons for calling a few preset configurations of find_multisyn using either a list of syri/aln files directly or a tsv containing this information
# For more information, see the find_multisyn docstring
def coresyn_from_tsv(path, **kwargs):
    return pasy.find_multisyn(*parse_input_tsv(path), only_core=True, **kwargs)
def crosssyn_from_tsv(path, **kwargs):
    return pasy.find_multisyn(*parse_input_tsv(path), only_core=False, **kwargs)
def coresyn_from_lists(syns, alns, **kwargs):
    return pasy.find_multisyn(syns, alns, only_core=True, **kwargs)
def crosssyn_from_lists(syns, alns, **kwargs):
    return pasy.find_multisyn(syns, alns, only_core=False, **kwargs)

def get_orgs_from_df(df):
    """Small utility function to get all organism from a DataFrame of `Pansyn` objects.
    :param df: A `DataFrame` containing `Pansyn` objects.
    :param_type df: `pandas.DataFrame`
    :returns: A `set` containing all of the organisms in a DataFrame of `pansyn` objects.
    """
    if df.empty:
        logger.error(f"get_orgs_from_df called with empty dataframe: {df}")
        raise ValueError("DF is empty!")
    return functools.reduce(lambda x, y: x.union(y), map(lambda x: set(x[1][0].ranges_dict.keys()), df.iterrows()))


def get_len(df):
    if df.empty:
        logger.error(f"get_len called with empty dataframe: {df}")
        raise ValueError("DF is empty!")
    return sum(map(lambda x: len(x.ref), map(lambda x: x[1][0], df.iterrows())))

# the warnings in the two functions below are spurious, see
# https://github.com/cython/cython/issues/1699
def tabularize_lens(df):
    if df.empty:
        logger.error(f"tabularize_lens called with empty dataframe: {df}")
        raise ValueError("DF is empty!")
    maxdegree = max(map(lambda x: x[1][0].get_degree(), df.iterrows()))
    return [sum(len(x[1][0].ref) for x in df.iterrows() if x[1][0].get_degree() == i+1) for i in range(maxdegree)]
    #return [sum(map(lambda x: len(x.ref), filter(lambda x: x.get_degree() == i + 1, map(lambda x: x[1][0], df.iterrows())))) for i in range(maxdegree)]

def tabularize_nos(df):
    if df.empty:
        logger.error(f"tabularize_nos called with empty dataframe: {df}")
        raise ValueError("DF is empty!")
    maxdegree = max(x[1][0].get_degree() for x in df.iterrows())
    return [sum(1 for x in df.iterrows() if x[1][0].get_degree() == i+1) for i in range(maxdegree)]

def get_stats(df):
    """Utility function to output some stats for a df containing computed pansyn objects.
    Calls get_len and tabularize_lens, prettyprints their output.
    """
    tot_len = get_len(df)
    if tot_len == 0:
        return "Empty!"

    lens = tabularize_lens(df)
    nos = tabularize_nos(df)
    avglens = list(map(lambda x: x[0]/x[1] if x[1] > 0 else 0, zip(lens, nos)))
    ret = f"Total syn length: {siprefix(tot_len)}\nDeg.\tTot. Length\tNo of Regions\tAvg. Length\n" + "\n".join([f"{i + 1}\t{siprefix(lens[i])}\t{nos[i]}\t{siprefix(avglens[i])}" for i, _ in enumerate(lens)])
    return ret

def siprefix(x: int):
    if x >= 1E9:
        return f"{x/1E9:.2f} Gbp"
    elif x >= 1E6:
        return f"{x/1E6:.2f} Mbp"
    elif x >= 1E3:
        return f"{x/1E3:.2f} kbp"
    else:
        return f"{x} bp"

def get_call_stats(syns, alns, **kwargs):
    """Utility function to call multisyn in a dataset and immediately compute the statistics using get_stats
    """
    df = pasy.find_multisyn(syns, alns, **kwargs)
    return get_stats(df)


def eval_combinations(syns, alns, cores=1):
    """Perform get_call_stats for all possible combinations of only_core and SYNAL
    """
    ret = ""
    for only_core, SYNAL in [(False, False), (False, True), (True, False), (True, True)]:
        ret += "core" if only_core else "all"
        ret += " pansynteny, "
        ret += "exact" if SYNAL else "approximate"
        ret += ":\n"
        ret += get_call_stats(syns, alns, only_core=only_core, SYNAL=SYNAL)
    return ret

def filter_pansyns(df, predicate):
    #TODO refactor to use Cpp vectors later
    # remember appropriate preallocation
    inds = df[0].apply(predicate)
    return df.loc[inds]

def apply_filtering(df, exp):
    return filter_pansyns(df, compile_filter(exp))

def filter_multisyn_df(df, rng, only_contained=False):
    """DEPRECATED, use filter_pansyns or apply_filtering instead
    Misc function for filtering a DF produced by find_multisyn for a certain range.
    Only the position on the reference is taken into account.
    Only the chromosome, start and end of the supplied `Range` are used, org and chromosome information is discarded.

    :param df: `find_multisyn` `DataFrame` of `Pansyn` objects.
    :type df: `DataFrame[Pansyn]`
    :param rng: `Range` for selecting the `Pansyn` objects.
    :param only_contained: switches between selecting any region intersecting or contained in the specified `Range`.
    :type only_contained: bool
    """
    def filter_fn(x):
        ref = x.ref
        # check if on same chr
        if not rng.chr == ref.chr:
            return False
        # check if contained:
        if rng.start < ref.start < rng.end and rng.start < ref.end < rng.end:
            return True
        # quit if only looking for contained regions
        if only_contained:
            return False
        # check if start or end within rng
        return rng.start < ref.start < rng.end or rng.start < ref.end < rng.end

    inds = df[0].apply(filter_fn)
    #print(inds)
    return df.loc[inds]

cpdef compile_filter_py(exp: str):
    return compile_filter(exp)

cdef compile_filter(exp: str):
    """
    Higher-Order function for compiling a filtering expression into a single, fast predicate.
    Returns a cdef lambda, check if this works properly when importing into python code, otherwise call in cpdef method.
    Runtime is in n*log(n) with n being len(exp) -- this could be made faster by using a proper parser, but I don't think this step is performance-limiting.

    # Ideas for DSL for filtering:
    # primitives: Range, degree >= number, len, chr, maybe alignment quality?
    # have preprocessor turn is pansyn into degree >= number of seqs
    # connections: and, or, xor, not
    """

    if len(exp) < 1:
        logger.warning("compile_filter called with empty string. This expression will match everything!")
        return lambda x: True

    # handle negation
    match = re.fullmatch("(\!|not)\s?(.*)", exp)
    if match:
        return lambda x: not compile_filter(match[2])(x)

    # handle and clauses
    match = re.fullmatch("\((.*)\)\s?(&|and)\s?\((.*)\)", exp, flags=re.IGNORECASE)
    if match:
        return lambda x: compile_filter(match[1])(x) and compile_filter(match[3])(x)

    # handle or clauses
    match = re.fullmatch("\((.*)\)\s?(\||or)\s?\((.*)\)", exp, flags=re.IGNORECASE)
    if match:
        return lambda x: compile_filter(match[1])(x) or compile_filter(match[3])(x)

    # handle xor clauses
    match = re.fullmatch("\((.*)\)\s?(\^|xor)\s?\((.*)\)", exp, flags=re.IGNORECASE)
    if match:
        return lambda x: compile_filter(match[1])(x) ^ compile_filter(match[3])(x)

    ## handle primitives
    # degree filters
    match = re.fullmatch("(deg|degree)\s?>=\s?(\d)+", exp, flags=re.IGNORECASE)
    if match:
        num = int(match[2])
        return lambda x: x.get_degree() >= num
    match = re.fullmatch("(deg|degree)\s?<=\s?(\d)+", exp, flags=re.IGNORECASE)
    if match:
        num = int(match[2])
        return lambda x: x.get_degree() <= num
    # maybe also implement > < etc

    # Len filters
    match = re.fullmatch("(len|length)\s?>=\s?(\d)+", exp, flags=re.IGNORECASE)
    if match:
        num = int(match[2])
        return lambda x: len(x.ref) >= num
    match = re.fullmatch("(len|length)\s?<=\s?(\d)+", exp, flags=re.IGNORECASE)
    if match:
        num = int(match[2])
        return lambda x: len(x.ref) <= num

    # organism filters
    match = re.fullmatch("(cont|contains)\s(.*)", exp, flags=re.IGNORECASE)
    if match:
        org = match[2]
        return lambda x: org in x.get_orgs()

    match = re.fullmatch("(contall|containsall)\s(.*)", exp, flags=re.IGNORECASE)
    if match:
        orgs = match[2].split(',')
        return lambda x: all([org.strip() in x.get_orgs() for org in orgs])

    match = re.fullmatch("(contany|containsany)\s(.*)", exp, flags=re.IGNORECASE)
    if match:
        orgs = match[2].split(',')
        return lambda x: any([org.strip() in x.get_orgs() for org in orgs])

    # handle position on reference, TODO maybe also do this for organism?
    match = re.fullmatch("(in)\s(.*)", exp, flags=re.IGNORECASE)
    if match:
        #TODO error handling?
        rng = Range.read_pff(None, match[2])
        return lambda x: x.ref in rng

    # chr filter
    match = re.fullmatch("(on)\s(.*)", exp, flags=re.IGNORECASE)
    if match:
        return lambda x: x.ref.chr == match[2]
    
    # find simple cases
    if exp.lower() == 'true':
        return lambda x: True
    elif exp.lower() == 'false':
        return lambda x: False
    # de-bracket expressions
    elif exp[0] == '(' and exp[-1] == ')': # no regex required!
        return compile_filter(exp[1:-1])



    logger.error(f"compile_filter called with invalid expression: {exp}")
    raise ValueError(f"compile_filter called with invalid expression: {exp}")

    

def filter_multisyn_df_chr(df, chr):
    """Misc function for filtering a DF produced by find_multisyn for a certain chromosome.
    Does essentially the same thing as `filter_multsyn_df`, but only uses chromosome information

    :param df: `find_multisyn` `DataFrame` of `Pansyn` objects.
    :type df: `DataFrame[Pansyn]`
    :param chr: Chromosome to select.
    :type chr: `str`
    """
    return df.loc[df[0].apply(lambda x: chr == x.ref.chr)]

def length_compare(syns, alns, cores=1):
    syns, alns = list(syns), list(alns)
    for _ in range(len(syns)):
        eval_combinations(syns, alns)
        syns = syns[1:]
        alns = alns[1:]

def pff_to_file(df, path):
    """Convenience wrapper for to_format to save to a file directly
    """
    with open(path, 'wt') as f:
        pasy.io.to_pff(df, f)

def pff_to_string(df, save_cigars=False):
    """Convenience wrapper for to_format, saves to a stringbuffer, returns the string.
    Mainly meant for printing small-ish callsets.
    """
    with io.StringIO() as buf:
        pasy.io.to_pff(df, buf, save_cigars=save_cigars)
        return buf.get_value()

