# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import numpy as np
import pandas as pd
from scipy.stats import *

from datetime import datetime
from multiprocessing import Pool
from functools import partial
from collections import deque, defaultdict, OrderedDict
from gzip import open as gzopen
from gzip import BadGzipFile

from collections import deque
import sys
import os
import logging
import psutil
import pysam
import re
from gc import collect

from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map
from libcpp.deque cimport deque as cpp_deq
cimport numpy as np
cimport cython

from pansyri.classes.coords import Range
from pansyri.classes.coords import Pansyn
from pansyri.classes.vars import SNV
import pansyri.util as util
import pansyri.classes.cigar as cigar

logger = util.CustomFormatter.getlogger(__name__)
np.random.seed(1)


### BEGIN func SECTION
## copied over from func

def cgtpl(cg):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    return [i.split(';') for i in cg.split(',')[:-1]]

# TODO maybe use https://pypi.org/project/pyfasta/ instead
def readfasta(f):
    out = {}
    chrid = ''
    chrseq = deque()

    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False

    try:
        if isgzip:
            with gzopen(f, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                            chrseq = deque()
                        else:
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        if chrid in out.keys():
                            logger.error(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                            raise ValueError()
                    else:
                        chrseq.append(line.strip().decode())
        else:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                        if chrid in out.keys():
                            logger.error(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                            raise ValueError()
                    else:
                        chrseq.append(line.strip())
    except Exception as e:
        raise Exception(e)

    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out

### END func SECTION


def samtocoords(f):
    logger = logging.getLogger('SAM reader')
    rc = {}        # Referece chromosomes
    rcs = {}        # Selected chromosomes
    al = deque()    # Individual alignment
    try:
        with open(f, 'r') as fin:
            for l in fin:
                if l[:3] == '@SQ':
                    c, s = 0, 0
                    for h in l.strip().split()[1:]:
                        h = h.split(':')
                        if h[0] == 'SN': c = h[1]
                        if h[0] == 'LN': s = int(h[1])
                    rcs[c] = s
                    continue
                elif l[0] == '@': continue

                l = l.split('\t')[:6]
                # if l[1] == '2064': break
                if l[2] == '*':
                    logger.warning(l[0]+ ' do not align with any reference sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')  # Skip rows corresponding to non-mapping sequences (contigs/scaffolds)
                    continue

                if 'M' in l[5]:
                    logger.error('Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: ' + l[5])
                    sys.exit()
                cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in l[5].replace('S', ';S,').replace('H', ';H,').replace('=', ';=,').replace('X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[:-1]]]
                if len(cgt) > 2:
                    if True in [True if i[1] in ['S', 'H'] else False for i in cgt[1:-1]]:
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                        sys.exit()

                bf = '{:012b}'.format(int(l[1]))

                rs = int(l[3])
                re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':    # forward alignment
                    if cgt[0][1] == '=':
                        qs = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qs = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qs = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qs = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qs, qe = qe, qs

                al.append([
                    rs,
                    re,
                    qs,
                    qe,
                    abs(re-rs) + 1,
                    abs(qs-qe) + 1,
                    format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                        [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                    1,
                    1 if bf[7] == '0' else -1,
                    l[2],
                    l[0],
                    "".join([str(i[0])+i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
                ])
                rcs[l[2]] = 1
            rcs = list(rcs.keys())
            for k in list(rc.keys()):
                if k not in rcs: logger.warning(l[0]+ ' do not align with any query sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')
    except Exception as e:
        logger.error('Error in reading SAM file: ' + str(e))
        sys.exit()
    al = DataFrame(list(al))
    al[6] = al[6].astype('float')
    al.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
    al.index = range(len(al.index))
    return al
# END

def readSAMBAM(fin, type='B'):
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        if type == 'B':
            findata = pysam.AlignmentFile(fin,'rb')
        elif type == 'S':
            return samtocoords(fin)
        else:
            raise ValueError("Wrong parameter")
    except ValueError as e:
        logger.error("Error in opening BAM/SAM file. " + str(e))
        sys.exit()
    except OSError as e:
        logger.error("Error in reading input file." + str(e))
        sys.exit()
    except Exception as e:
        logger.error("Unexpected error in opening BAM/SAM file. " + str(e))
        sys.exit()

    try:
        qry_prim = {}
        ref_prim = {}
        cgdict = {1:'I', 2:'D', 7:'=', 8:'X'}
        coords = {}
        index = 0
        for aln in findata:
            index += 1
            ## Check whether every sequence has at least one primary alignment
            if aln.reference_name is not None:
                if aln.reference_name not in ref_prim.keys():
                    ref_prim[aln.reference_name] = False
            if aln.query_name not in qry_prim.keys():
                qry_prim[aln.query_name] = False
            if aln.reference_name is not None:
                if not ref_prim[aln.reference_name]:
                    if aln.flag < 256:
                        ref_prim[aln.reference_name] = True
            if not qry_prim[aln.query_name]:
                if aln.flag < 256:
                    qry_prim[aln.query_name] = True

            ## Pass non-alinging chromosomes
            if aln.cigarstring is None:
                logger.warning(aln.query_name + ' do not align with any reference chromosome and cannot be analysed')
                continue

            ## Check CIGAR:
            if False in [False if i[0] not in [1,2,4,5,7,8] else True for i in aln.cigartuples]:
                logger.error("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(aln.cigarstring))
                sys.exit()
            if len(aln.cigartuples) > 2:
                if True in [True if i[0] in [4,5] else False for i in aln.cigartuples[1:-1]]:
                    logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                    sys.exit()

            ## Parse information from the aln object
            astart = aln.reference_start+1
            aend = aln.reference_end
            is_inv = True if np.binary_repr(aln.flag,12)[7] == '1' else False
            if not is_inv:
                if aln.cigartuples[0][0] in [4,5]:
                    bstart = aln.cigartuples[0][1]+1
                else:
                    bstart = 1
                bend = bstart + aln.query_alignment_length - 1
            else:
                if aln.cigartuples[-1][0] in [4,5]:
                    bend = aln.cigartuples[-1][1]+1
                else:
                    bend = 1
                bstart = bend + aln.query_alignment_length - 1
            alen = abs(aend - astart) + 1
            blen = abs(bend - bstart) + 1
            iden = format((sum([i[1] for i in aln.cigartuples if i[0] == 7])/sum([i[1] for i in aln.cigartuples if i[0] in [1,2,7,8]]))*100, '.2f')
            adir = 1
            bdir = -1 if is_inv else 1
            achr = aln.reference_name
            bchr = aln.query_name
            seq = aln.query_sequence
            cg = "".join([str(i[1]) + cgdict[i[0]] for i in aln.cigartuples if i[0] not in [4,5]])
            coords[index] = [astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg, seq]

        ## Give warning for chromosomes which do not have any primary alignment
        for k,v in ref_prim.items():
            if not v:
                logger.warning('No primary alignment found for reference sequence ' + k +'. This could mean that the entire chromosome '+ k +' is repeated.')
        for k,v in qry_prim.items():
            if not v:
                logger.warning('No primary alignment found for query sequence ' + k +'. This could mean that the entire chromosome '+ k + ' is repeated.')

        ## Return alignments
        coords = pd.DataFrame.from_dict(coords, orient= 'index')
        coords.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        coords.columns = ["astart", "aend", "bstart", "bend", "alen", "blen", "iden", "adir", "bdir", "achr", "bchr", "cg", "seq"]
        return coords
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        sys.exit()
# END

def readPAF(paf):
    coords = deque()
    logger = logging.getLogger('Reading PAF file')
    try:
        with open(paf, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                astart = int(line[7]) + 1
                aend = int(line[8])
                adir = 1
                bdir = 1 if line[4] == '+' else -1
                bstart = int(line[2]) + 1 if bdir == 1 else int(line[3])
                bend = int(line[3]) if bdir == 1 else int(line[2]) + 1
                alen = abs(aend - astart) + 1
                blen = abs(bend - bstart) + 1 if bdir == 1 else bstart - bend + 1
                cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg']
                if len(cg) != 1:
                    logger.error("CIGAR string is not present in PAF at line {}. Exiting.".format("\t".join(line)))
                    sys.exit()
                cg = cg[0]
                ## Check CIGAR:
                if not all([True if i[1] in {'I', 'D', 'H', 'S', 'X', '='} else False for i in cgtpl(cg)]):
                    logger.error("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(cg))
                    sys.exit()
                if len(cgtpl(cg)) > 2:
                    if any([True if i[1] in {'H', 'S'} else False for i in cgtpl(cg)]):
                        logger.error("Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + str(cg))
                        sys.exit()

                iden = round((sum([int(i[0]) for i in cgtpl(cg) if i[1] == '='])/sum([int(i[0]) for i in cgtpl(cg) if i[1] in {'=', 'X', 'D', 'I'}]))*100, 2)
                achr = line[5]
                bchr = line[0]
                coords.append([astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg])
        coords = pd.DataFrame(coords)
        coords.sort_values([9,0,1,2,3,10], inplace = True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        coords.columns = ["astart", "aend", "bstart", "bend", "alen", "blen", "iden", "adir", "bdir", "achr", "bchr", "cg"]
        return coords
    except FileNotFoundError:
        logger.error("Cannot open {} file. Exiting".format(paf))
        sys.exit()
    except ValueError as e:
        logger.error("Error in reading PAF: {}. Exiting".format(e))
        sys.exit()
# END

# pasted from plotsr, parsing syri output
VARS = ['SYN', 'SYNAL', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
cpdef readsyriout(f):
    # Reads syri.out. Select: achr, astart, aend, bchr, bstart, bend, srtype
    logger = logging.getLogger("readsyriout")
    syri_regs = deque()
    skipvartype = ['CPG', 'CPL', 'DEL', 'DUPAL', 'HDR', 'INS', 'INVAL', 'INVDPAL', 'INVTRAL', 'NOTAL', 'SNP', 'TDM', 'TRANSAL']
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            # TODO: DECIDE WHETHER TO HAVE STATIC VARS OR FLEXIBLE ANNOTATION
            if l[10] in VARS:
                syri_regs.append(l)
            else:
                if l[10] not in skipvartype:
                    skipvartype.append(l[10])
                    logger.warning("{} is not a valid annotation for alignments in file {}. Alignments should belong to the following classes {}. Skipping alignment.".format(l[10], f, VARS))

    try:
        df = pd.DataFrame(list(syri_regs))[[0, 1, 2, 5, 6, 7, 10]]
    except KeyError:
        raise ImportError("Incomplete input file {}, syri.out file should have 11 columns.".format(f))
    df[[0, 5, 10]] = df[[0, 5, 10]].astype(str)
    try:
        df[[1, 2, 6, 7]] = df[[1, 2, 6, 7]].astype(int)
    except ValueError:
        raise ValueError("Non-numerical values used as genome coordinates in {}. Exiting".format(f))
    # chr ID map
    chrid = []
    chrid_dict = OrderedDict()
    for i in np.unique(df[0]):
        chrid.append((i, np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]))
        chrid_dict[i] = np.unique(df.loc[(df[0] == i) & (df[10] == 'SYN'), 5])[0]
    df.columns = ['achr', 'astart', 'aend', 'bchr', 'bstart', 'bend',  'type']
    return df, chrid_dict

cpdef extract_syri_snvs(fin):
    syri_regs = deque()
    with open(f, 'r') as fin:
        for line in fin:
            l = line.strip().split()
            if l[10] == 'SNP':
                #TODO maybe store annotation information from fields 8-10
                snv = SNV(Position('a', 'x', l[0], int(l[1])), Position('b', 'x', l[5], int(l[6])), l[4], l[5])
                syri_regs.append(SNV)

    df = pd.DataFrame(list(syri_regs))#[[0, 1, 3, 4, 5, 6, 8, 9, 10]]
    #TODO maybe do chromosome mapping?
    return df

HEADER="""##INFO=<ID=END,Number=1,Type=Integer,Description="End position on reference genome">
##ALT<ID=CORESYN,Description="Core syntenic region (syntenic between any two samples)">
##ALT<ID=CROSSSYN,Description="Cross syntenic region (syntenic between any two samples for a strict subset of the samples)>
##FORMAT=<ID=CHR,Number=1,Type=String,Description="Chromosome in this sample">
##FORMAT=<ID=START,Number=1,Type=Integer,Description="Start position in this sample">
##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position  in this sample">
##FORMAT=<ID=SYN,Number=1,Type=Integer,Description="1 if this region is syntenic to reference, else 0">"""
##FORMAT=<ID=HAP,Number=1,Type=Character,Description="Unique haplotype identifier">"""

cpdef extract_syntenic_from_vcf(syns, inpath, outpath, force_index=True, org='ref', ref=None, keep_nonsyn_calls=False):
    """
    Extract syntenic annotations from a given VCF.
    A tabix-indexed VCF is required for this; by default, the input VCF is reindexed (and gzipped) with the call.
    If the supplied VCF already has a tabix index, `force_index` may be set to false.
    """
    cdef:
        vcfin = pysam.VariantFile(inpath)
        vcfout = pysam.VariantFile(outpath, 'w', header=vcfin.header)
        orgs = util.get_orgs_from_df(syns)
        header_chrs = set(vcfin.header.contigs)
        int crosscounter = 0
        int corecounter = 0

    if not set(vcfin.header.samples).issubset(orgs):
        logger.warning("Input VCF contains organisms not in PFF file! Double-Check names used in .tsv. Truncating VCF.")
        vcfin.subset_samples(orgs)

    # read reference if it hasn't been read already
    if type(ref) != dict:
        logger.info("Reading in Reference Fasta")
        ref = readfasta(ref)


    orgsvcf = list(vcfin.header.samples) # select only 

    # force indexing to allow for calling fetch later.
    if force_index and not vcfin.index_filename:
        vcfin.close()
        pysam.tabix_index(inpath, force=True, preset='vcf', keep_original=True)
        inpath += ".gz" # no way to turn off automatic compression, apparently
        vcfin = pysam.VariantFile(inpath)

    # add header required for storing PANSYN annotations
    for line in HEADER.splitlines():
        vcfout.header.add_line(line)

    # add pansyn regions and all variation therein
    for syn in syns.iterrows():
        syn = syn[1][0]
        rng = syn.ref if org == 'ref' else syn.rngs[org]
        rec = vcfout.new_record()
        rec.start = rng.start
        rec.pos = rec.start
        rec.stop = rng.end
        chrom = rng.chr

        if chrom not in header_chrs:
            if ref:
                # add length if it is known from the reference
                vcfout.header.add_line("##contig=<ID={},length={}>".format(chrom, len(ref[chrom])))
            else:
                vcfout.header.add_line("##contig=<ID={}>".format(chrom))

        rec.chrom = chrom
        if syn.get_degree() == len(orgs):
            if ref:
                rec.alleles = [ref[rec.chrom][rec.start], "<CORESYN>"]
            else:
                rec.alleles = ["<SYN>", "<CORESYN>"]
            rec.id = "CORESYN{}".format(corecounter)
            corecounter += 1
        else:
            if ref:
                rec.alleles = [ref[rec.chrom][rec.start], "<CROSSSYN>"]
            else:
                rec.alleles = ["<SYN>", "<CROSSSYN>"]
            rec.id = "CROSSSYN{}".format(crosscounter)
            crosscounter += 1

        # write the pansyn annotation
        for org in orgsvcf:
            if org in syn.get_orgs():
                rng = syn.ranges_dict[org]
                rec.samples[org].update({'SYN':1, 'CHR':rng.chr, 'START': rng.start, 'END': rng.end})
            else:
                rec.samples[org].update({'SYN': 0})
        vcfout.write(rec)

        # write the small variants in the pansyn region
        for rec in vcfin.fetch(rng.chr, rng.start, rng.end + 1): # pysam is half-inclusive
            # iterate through organisms, remove any data that is not syntenic
            if not keep_nonsyn_calls:
                for org in orgsvcf:
                    if org not in syn.get_orgs():
                        del rec[org]
            vcfout.write(rec) # this is failing, but still writing the correct output? WTF?

    #vcfout.close()
    #vcfin.close()



cpdef extract_syri_regions(fin, ref='a', anns=['SYN'], reforg='ref', qryorg='qry'):
    """
    Given a syri output file, extract all regions matching a given annotation.
    """
    # columns to look for as start/end positions
    refchr = ref + "chr"
    refhaplo = "x"
    refstart = ref + "start"
    refend = ref + "end"

    qry = 'b' if ref == 'a' else 'a' # these seem to be the only two values in syri output
    qrychr = qry + "chr"
    qryhaplo = "x"
    qrystart = qry + "start"
    qryend = qry + "end"

    buf = deque()
    raw, chr_mapping = readsyriout(fin) #TODO? handle chr_mapping
    raw = pd.concat([raw.loc[raw['type'] == ann] for ann in anns])
    # if implementing filtering later, filter here

    for row in raw.iterrows():
        row = row[1]
        buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
            Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
            ])

    return pd.DataFrame(data=list(buf), columns=[reforg, qryorg])

def extract_syri_regions_to_list(fins, qrynames, cores=1, **kwargs):
    """
    `extract_syri_regions`, but for processing a list of inputs
    """
    if len(fins) != len(qrynames):
        logger.error(f"Infiles and qrynames lists lengths not matching. Offending lists: {fins} and {qrynames}")
    partial = lambda x, qryname: extract_syri_regions(x, qryorg=qryname, **kwargs)

    if cores == 1:
        syns = [partial(fin, qryname) for fin, qryname in zip(fins, qrynames)]
    else:
        with Pool(cores) as pool:
            syns = pool.map(partial, fins)

    return syns
    #return [extract_syri_regions(fin, **kwargs,\
    #        #reforg=fin.split('/')[-1].split('_')[0],\
    #        qryorg=fin.split('/')[-1].split('_')[-1].split('syri')[0])\
    #        for fin in fins]

cpdef save_to_vcf(syns, outf, ref=None, cores=1):
    #TODO add functionality to incorporate reference information as optional argument
    cdef:
        out = pysam.VariantFile(outf, 'w')
        int corecounter = 1 # 1-based region indexing
        int crosscounter = 1
        # ensure consistent, alphabetical sorting of organisms
        orgs = sorted(util.get_orgs_from_df(syns))
        int orgsc = len(orgs)
        header_chrs = set() # do dynamically in python, hopefully more efficiently than looping twice

    # prepare appropriate header file
    for line in HEADER.splitlines():
        out.header.add_line(line)

    if type(ref) != dict:
        logger.info("Reading in Reference Fasta")
        ref = readfasta(ref)

    #out.header.add_samples(util.get_orgs_from_df(syns)) # according to the documentation, this works, but the function doesn't seem to exist...
    for org in orgs:
        out.header.add_sample(org)

    # add each pansyn object
    for syn in syns.iterrows():
        syn = syn[1][0]

        rec = out.new_record()
        # instantiate empty, then fill later
        # instantiating with keyword arguments is unstable according to the documentation
        rec.start = syn.ref.start
        rec.pos = syn.ref.start
        ## Chr needs to be a number, format it:
        #match = re.fullmatch(r"\D*?(\d+)\D*", syn.ref.chr)
        #if not match:
        #    logger.error("VCF exporting only accepts chr names only containing one number such as Chr12, but not chr names containing more than one number, e.g. Chr12_1! Offending chr name:" + syn.ref.chr)
        #else:
        #    chrom = match[1]

        ## store Chr as string for now, maybe change later
        chrom = syn.ref.chr
        if chrom not in header_chrs:
            if ref:
                # add length if it is known from the reference
                out.header.add_line("##contig=<ID={},length={}>".format(chrom, len(ref[chrom])))
            else:
                out.header.add_line("##contig=<ID={}>".format(chrom))

        rec.chrom = chrom

        rec.stop = syn.ref.end # apparently this exists? what does it do?
        if syn.get_degree() == orgsc:
            if ref:
                rec.alleles = [ref[rec.chrom][rec.start], "<CORESYN>"]
            else:
                rec.alleles = ["<SYN>", "<CORESYN>"]
            rec.id = "CORESYN{}".format(corecounter)
            corecounter += 1
        else:
            if ref:
                rec.alleles = [ref[rec.chrom][rec.start], "<CROSSSYN>"]
            else:
                rec.alleles = ["<SYN>", "<CROSSSYN>"]
            rec.id = "CROSSSYN{}".format(crosscounter)
            crosscounter += 1

        # input the values for every organism
        for org in orgs:
            if org in syn.get_orgs():
                rng = syn.ranges_dict[org]
                ## comment out chr to int conversion for now
                # Chr needs to be a number, format it:
                #match = re.fullmatch(r"\D*?(\d+)\D*", syn.ref.chr)
                #chrom = 1
                #if not match:
                #    logger.error("VCF exporting only accepts chr names only containing one number such as Chr12, but not chr names containing more than one number, e.g. Chr12_1! Offending chr name:" + syn.ref.chr)
                #    rec.samples[org].update({'SYN':1, 'START': rng.start, 'END': rng.end})
                #    continue
                #else:
                #    chrom = int(match[1])
                rec.samples[org].update({'SYN':1, 'CHR':rng.chr, 'START': rng.start, 'END': rng.end})
            else:
                rec.samples[org].update({'SYN': 0})

        out.write(rec)
    out.close()

cpdef save_to_pff(df, buf, save_cigars=True):
    """Takes a df containing `Pansyn` objects and writes them in pansyri file format to `buf`.
    Can be used to print directly to a file, or to print or further process the output.
    """
    # output organisms in lexicalic ordering
    orgs = sorted(util.get_orgs_from_df(df))
    buf.write("#ANN\tref\t")
    buf.write("\t".join(orgs))

    for row in df.iterrows():
        pansyn = row[1][0]
        buf.write("\nSYN\t") # only handle SYNs for now
        buf.write(pansyn.ref.to_pff())
        for org in orgs:
            buf.write("\t")
            if org in pansyn.ranges_dict:
                buf.write(pansyn.ranges_dict[org].to_pff())
                if save_cigars and pansyn.cigars_dict:
                    buf.write(",")
                    buf.write(pansyn.cigars_dict[org].to_string())
            else:
                buf.write(".")

    buf.write("\n")

cpdef read_pff(f):
    """Takes a file object or path to a file in PFF format and reads it in as a DataFrame.
    """
    syns = deque()
    orgs = f.readline().strip()[1:].split("\t")[2:] # 0 is ANN, 1 is ref
    for l in f:
        l = l.strip().split('\t')
        if l[0] == 'SYN': # line contains a pansyn region
            syn = Pansyn(Range.read_pff("ref", l[1]), # extract reference range
                {orgs[i]:Range.read_pff(orgs[i], cell.split(",")[0]) # extract ranges dict
                        for i, cell in enumerate(l[2:]) if cell != '.'},
                {orgs[i]:cigar.cigar_from_string(cell.split(',')[1])
                 for i, cell in enumerate(l[2:]) if cell != '.' and len(cell.split(',')) > 1} # extract cigars dict
            )
            syns.append(syn)


    return pd.DataFrame(data=list(syns)) # shouldn't require sorting
