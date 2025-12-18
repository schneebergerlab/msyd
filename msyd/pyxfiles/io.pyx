# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import numpy as np
import pandas as pd
from scipy.stats import *
import pysam
import functools

#from multiprocessing import Pool
from collections import deque, defaultdict, OrderedDict
from gzip import open as gzopen
from gzip import BadGzipFile

from collections import deque
import sys
import os
import logging
#import re

cimport numpy as np

from msyd.syngraph import Node, trace_org
from msyd.coords import Range, read_psf_range, Panco
from msyd.multisyn import Multisyn, MultisynContainer, ChromContainer
from msyd.orgs import OrgContainer
from msyd.vars import SNV
import msyd.util as util
import msyd.cigar as cigar

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

                rstart = int(l[3])
                rend = rstart - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':    # forward alignment
                    if cgt[0][1] == '=':
                        qstart = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qstart = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qend = qstart - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qstart = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qstart = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qend = qstart - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qstart, qend = qend, qstart

                al.append([
                    rstart,
                    rend,
                    qstart,
                    qend,
                    abs(rend-rstart) + 1,
                    abs(qstart-qend) + 1,
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
                raise ValueError("Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(aln.cigarstring))
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

alnfilelookup = {
        'sam': readSAMBAM,
        'bam': readSAMBAM,
        'paf': readPAF
        }

def split_alndf_by_chrom(alndf, chromid="achr"):
    """
    Takes a DF of alignments, returns a Dictionary mapping chromosome names to alignments on them.
    As Chromosome names, the contents of the `chromid` arg are taken.
    Fairly inefficient, would be faster to do this already while reading in the alns.
    """
    return {chrom: df for chrom, df in alndf.groupby(by=chromid)}

def collate_by_chrom(alndfs, chromid="achr"):
    out = defaultdict(list)
    for alndf in alndfs:
        for chrom, alns in split_alndf_by_chrom(alndf, chromid=chromid).items():
            out[chrom].append(alns)
    return out

cpdef read_alnsfile(fin):
    """
    Reads in pairwise all vs all alignments.
    The file containing each alignment is in the third column, the first two columns contain the samples used as reference and alternative in the alignment.
    The fourth column may optionally contain the filetype (sam, bam or paf).
    If no fourth column is present, the filetype will be inferred from the file ending
    Empty lines and lines starting with '#' are ignored
    """
    if isinstance(fin, str):
        fin = open(fin, 'rt')

    out = dict()
    for line in fin:
        line = line.strip()
        # ignore empty lines and comments
        if line == '' or line[0] == '#':
            continue
        
        # parse TSV
        line = line.split('\t')
        ref = line[0].strip()
        alt = line[1].strip()
        path = line[2].strip()

        # read in aln file
        ftype = line[3].strip() if len(line) > 3 else path.split('.')[-1]
        aln = alnfilelookup[ftype.lower()](path)
        
        # add to dict of dicts structure
        if ref not in out:
            out[ref] = {alt: aln}
        else:
            if alt in out[ref]: # warn on duplicate
                logger.warning(f"Duplicate alignment from {alt} to {ref}. Using last one!")
            out[ref][alt] = aln

    # check if all pairwise alignments are present, otherwise warn
    for org in out:
        others = set(out)
        others.remove(org)
        if not set(out[org]) == others:
            logger.warning(f"Not all pairwise alignments present for {org}! This may cause errors during realignment.")

    return out


cpdef void save_to_vcf(chromcont: ChromContainer, outf: Union[str, os.PathLike], ref=None, cores=1, add_cigar=False, add_identity=True):
    #TODO add functionality to incorporate reference information as optional argument
    cdef:
        out = pysam.VariantFile(outf, 'w')
        counter = Panco.start_counter("") # add chrom later
        # ensure consistent, alphabetical sorting of organisms
        int orgsc = len(orgs)
        header_chrs = set() # do dynamically in python, hopefully more efficiently than looping twice

    # prepare appropriate header file
    for line in HEADER.splitlines():
        out.header.add_line(line)

    if ref and type(ref) != dict:
        logger.info("Reading in Reference Fasta")
        ref = readfasta(ref)
    elif not ref:
        logger.warning("No Reference specified, not saving Ref Sequence in VCF!")

    #out.header.add_samples(util.get_orgs_from_df(syns)) # according to the documentation, this works, but the function doesn't seem to exist...
    for org in iter(orgs):
        out.header.add_sample(org)

    # add each multisyn object
    for syn in iter(chromcont):

        rec = out.new_record()
        # instantiate empty, then fill later
        # instantiating with keyword arguments is unstable according to the documentation
        rec.start = syn.ref.start
        rec.pos = syn.ref.start

        ## store Chr as string for now, maybe change later
        chrom = syn.ref.chr
        counter.chrom = chrom # make sure to update it
        if chrom not in header_chrs:
            #logger.info(f"save_to_vcf Adding {chrom} to header")
            header_chrs.add(chrom)
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
            rec.id = "CORESYN:{}".format(counter)
            counter.increment_c()
        else:
            if ref:
                rec.alleles = [ref[rec.chrom][rec.start], "<MERASYN>"]
            else:
                rec.alleles = ["<SYN>", "<MERASYN>"]
            rec.id = "MERASYN:{}".format(counter)
            counter.increment_m()

        #rec.info['NS'] = syn.get_degree() # update NS column, include not only orgs in sample now

        # input the values for every organism
        for org in iter(orgs):
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

                if syn.cigars_dict:
                    cg = syn.cigars_dict[org]
                    if add_cigar:
                        rec.samples[org].update({'CG': cg.to_string()})
                    if add_identity:
                        rec.samples[org].update({'AI': int(cg.get_identity()*100)})
            else:
                rec.samples[org].update({'SYN': 0})

        out.write(rec)
    out.close()

cpdef save_to_psf(chromcont, buf, save_cigars=True, force_ref_pos=False, ref="ref"):
    """
    Takes a `ChromContainer` object containing one `MultisynContainer` per chromosome and writes them to buf.
    Preserves the sorting of the DFs, sorts chroms lexicallicaly.
    Calls to `save_cont_to_psf`.
    """
    if chromcont.is_empty():
        raise ValueError("Empty dfmap provided!")

    chromcont.orgs.sort()

    # write header; assumes the first chrom contains all orgs at least once
    buf.write("#CHR\tSTART\tEND\tANN\tREP\tRCHR\tRSTART\tREND\t")
    buf.write("\t".join(chromcont.orgs.get_names()))
    buf.write("\n")

    ## write records
    _save_msyncont_to_psf = functools.partial(save_msyncont_to_psf, buf=buf, emit_header=False, save_cigars=save_cigars, force_ref_pos=force_ref_pos, ref=ref)
    chromcont.apply_chroms(_save_msyncont_to_psf)
    #TODO parallelize?
    #TODO print comment about which chrom is starting?
    #save_msyncont_to_psf(iter(chromcont), buf, chromcont.orgs, emit_header=False, save_cigars=save_cigars, force_ref_pos=force_ref_pos)

cpdef save_msyncont_to_psf(chrom, msyncont, buf, save_cigars=True, emit_header=True, force_ref_pos=False, ref="ref"):
    """Takes a  a `MultisynContainer` per chromosome and writes them in population synteny file format to `buf`.
    Can be used to print directly to a file, or to print or further process the output.
    """
    cdef:
        #int n = len(msyncont.orgs)
        str coreend = "0"
        object counter = Panco.start_counter(chrom)

    if emit_header:
        buf.write("#CHR\tSTART\tEND\tANN\tREP\tRCHR\tRSTART\tREND\t")
        buf.write("\t".join(iter(msyncont.orgs)))
        buf.write("\n")

    for core, mesyns in msyncont.iter_cores_acc():
        # first, write non-ref-position merasynteny
        # write to the first position it can be
        # maybe this should be annotated for the entire range it can be instead (coreend+1:syn.start-1)
        for mesyn in mesyns:
            counter.increment_m()
            # write the BED-like pre record cols
            panco_id = f"MERASYN{counter}" if mesyn.get_degree() > 1 else f"PRIVATE{counter}"

            if mesyn.ref.org == ref:
                buf.write('\t'.join([mesyn.ref.chr, str(mesyn.ref.start), str(mesyn.ref.end), panco_id, mesyn.ref.org, '.', '.', '.', '']))
            elif force_ref_pos:
                buf.write('\t'.join([chrom, coreend, coreend, panco_id, mesyn.ref.org, mesyn.ref.chr, str(mesyn.ref.start), str(mesyn.ref.end), '']))
            else:
                buf.write('\t'.join(['.', '.', '.', panco_id, mesyn.ref.org, mesyn.ref.chr, str(mesyn.ref.start), str(mesyn.ref.end), '']))
            # write the record
            write_multisyn(mesyn, buf, msyncont.orgs, save_cigars=save_cigars)

        # write coresyn region
        if core:
            counter.increment_c()
            buf.write('\t'.join([core.ref.chr, str(core.ref.start), str(core.ref.end), f"CORESYN{counter}", core.ref.org, '.', '.', '.', '']))
            write_multisyn(core, buf, msyncont.orgs, save_cigars=save_cigars)
            coreend = str(core.ref.end + 1)

    buf.write("\n")
    buf.flush()

cdef write_multisyn(multisyn, buf, orgs, save_cigars=False):
    """Function to write a multisyn object in PSF style to buf.
    Does not write the BED-like first part of the annotation.
    :param multisyn: multisyn object to write
    :param buf: buffer to write to
    :param orgs: ordering of organisms to use (should be sorted)
    """
    buf.write('\t'.join([(multisyn.ranges_dict[org].to_psf()\
                if not (save_cigars and multisyn.cigars_dict) else\
                    ','.join([multisyn.ranges_dict[org].to_psf(), multisyn.cigars_dict[org].to_string()]) )
                         if (not multisyn.is_private()) and (org in multisyn.ranges_dict) else
                         (multisyn.ref.to_psf() if multisyn.ref.org == org else '.') # if there is no synteny, put a .
                 for org in iter(orgs)])
              )
    buf.write("\n")

cpdef object read_psf(fin): # -> ChromContainer
    """
    Takes a file object or path to a file in PSF format and reads it in as a DataFrame of Multisynteny objects.
    Supports the new version of PSF format; for legacy files, use the deprecated version of this function.
    """
    if isinstance(fin, str):
        fin = open(fin, 'rt')

    #CHR  START  END  ANN  REF  CHR  START  END  G1  G2  G3...
    cdef:
        line = fin.readline().strip().split() # Header line
        orgs = OrgContainer.from_list(line[8:]) # preserve file order
        chromdict = dict()#defaultdict(lambda: MultisynContainer(orgs))

    for line in fin:
        line = line.strip().split()
        if line == []: continue

        chrom = line[0]
        reforg = line[4]

        refrng = Range('ref', line[0], int(line[1]), int(line[2])) if reforg == 'ref'\
            else Range(reforg, line[5], int(line[6]), int(line[7]))

        syn = Multisyn(refrng, {}, None)

        for org, entry in zip(iter(orgs), line[8:]):
            if entry == '.' or org == reforg: # skip empty records and ref
                continue

            vals = entry.split(',')
            syn.ranges_dict[org] = read_psf_range(vals[0], org=org)
            
            # read cigars if present
            if len(vals) > 1:
                if syn.cigars_dict:
                    syn.cigars_dict[org] = cigar.cigar_from_string(vals[1])
                else: # initialise if it hasn't been already
                    syn.cigars_dict = {org: cigar.cigar_from_string(vals[1])}
        # add read in syn to output
        if chrom not in chromdict:
            chromdict[chrom] = MultisynContainer(orgs)
        chromdict[chrom].append(syn)

    # clean up, return
    fin.close()

    # return as ChromContainer
    return ChromContainer(chromdict, orgs)

cpdef save_to_gfa1(chromcont, buf, rgfa_tags=True, vg_header=True, tag_orgs_s=True, tag_orgs_l=True, walks_orgs=True):
    """
    Takes a map of chrom IDs to DFs containing Node objects and writes them to buf in GFA1 format.
    Preserves the sorting of the DFs (should be in topological sorting) and sorts chroms lexicallicaly.
    Calls to `save_df_to_gfa1`.
    :params:
    :tag_orgs_s, tag_orgs_l, walks_orgs: True (default), False or Set of organism names. Whether to write all, none, or only the specified organisms as tags on the segments, links, or as walks at the end of the GFA.
    :rgfa_tags: True (default), False. Whether to emit the rGFA tags of `SN`, `SO` and `SR`. The latter will always be 1 by convention.
    :vg_header: True (default), False. Whether to write vg's extended headers to the output GFA1.
    """
    if chromcont.is_empty():
        raise ValueError("Empty dfmap provided!")

    # get list of all orgs to simplify regularization
    # start node contains all orgs
    logger.info(f"Orgs found: {chromcont.orgs}")
    # regularize input args to sets
    if tag_orgs_s == True or tag_orgs_s == "True": # == required to not match truthy nonempty sets
        tag_orgs_s = chromcont.orgs
        #TODO handle organism regularization, or just handle in main?
    elif tag_orgs_s == False or tag_orgs_s == "False":
        tag_orgs_s = OrgContainer()
    else:
        tag_orgs_s = OrgContainer.from_list(tag_orgs_s)


    if tag_orgs_l == True or tag_orgs_l == "True": # == required to not match truthy nonempty sets
        tag_orgs_l = chromcont.orgs
    elif tag_orgs_l == False or tag_orgs_l == "False":
        tag_orgs_l = OrgContainer()
    else:
        tag_orgs_l = OrgContainer.from_list(tag_orgs_l)

    if walks_orgs == True or walks_orgs == "True": # == required to not match truthy nonempty sets
        walks_orgs = chromcont.orgs
    elif walks_orgs == False or walks_orgs == "False":
        walks_orgs = OrgContainer()
    else:
        walks_orgs = OrgContainer.from_list(walks_orgs)

    logger.info(f"Tracing on S lines: {tag_orgs_s}")
    logger.info(f"Tracing on L lines: {tag_orgs_l}")
    logger.info(f"Computing W lines for {walks_orgs}")

    ## write header
    buf.write("H\tVN:Z:1.2")
    if vg_header and (tag_orgs_s or tag_orgs_l):
        buf.write("\tRS:Z:" + " ".join(tag_orgs_s + tag_orgs_l))
    buf.write("\n")

    ## write contents
    #NOTE could parallelize?
    _save_msyncont_to_gfa1 = functools.partial(save_msyncont_to_gfa1, buf=buf, rgfa_tags=True, tag_orgs_s=tag_orgs_s, tag_orgs_l=tag_orgs_l, walks_orgs=walks_orgs)
    chromcont.apply_chroms(_save_msyncont_to_gfa1)

    logger.info(f"Finished writing GFA")

cpdef save_msyncont_to_gfa1(chrom, msyncont, buf, tag_orgs_s=set(), tag_orgs_l=set(), rgfa_tags=True, walks_orgs=set()):
    # get start and end node from the beginning of the DF
    nodeiter = iter(msyncont)
    startnode = next(nodeiter)
    endnode = next(nodeiter)
    buf.write(f"# <chrom:{chrom}>\n")

    # write S and L lines corresponding to nodes
    for node in nodeiter:
        buf.write(node.to_gfa1(rgfa_tags=rgfa_tags, tag_orgs_s=tag_orgs_s, tag_orgs_l=tag_orgs_l))
        buf.write("\n")
        
    # write W lines
    if walks_orgs:
        logger.info(f"Writing Walks for {walks_orgs}")
        for org in walks_orgs:
            # example walk line from the GFA1 spec
            # W	NA12878	1	chr1	0	11	>s11<s12>s13
            # offload path tracing to fn in msyd.syngraph
            walk = trace_org(startnode, org)
            if not walk:
                logger.warning(f"Empty trace found for {org}!")
                continue
            startrng = walk[0].msyn.ranges_dict[org] if org in walk[0].msyn.ranges_dict else walk[0].msyn.ref
            # write fixed part of line
            buf.write(f"\nW\t{org}\t{startrng.start}\t{startrng.chr}")
            #TODO finish writing non-fixed part of line
            # notes
            # think if it makes sense to combine this with msyn refactor
            # => does including the ref in ranges_dict break stuff in intersection/realignment?

    buf.write(f"# </chrom:{chrom}>\n")
    logger.info(f"Finished {chrom} part of GFA")

cpdef read_old_psf(fin):
    """
    DEPRECATED, for reading PSF files produced by v0.2
    Takes a file object or path to a file in PSF format and reads it in as a DataFrame of Multisynteny objects.
    Supports the new version of PSF format; for legacy files, use the deprecated version of this function.
    """
    syns = deque()
    if isinstance(fin, str):
        fin = open(fin, 'rt')

    line = fin.readline().strip().split()
    samples = line[4:]
    for line in fin:
        # if line == '': continue
        # if line is None: continue
        line = line.strip().split()
        if line == []: continue
        #try:
        #    anno = line[3]
        #except IndexError:
        #    logger.error(f"Invalid line encountered while reading PSF: {line}")

        refrng = Range('ref', line[0], int(line[1]), int(line[2]))

        # split once to reuse in multisyn construction loop
        samplecells = [cell.split(';') for cell in line[4:]]

        # a single line may contain multiple merasyn records if the PSF is collapsed
        for i in range(len(samplecells[0])): # will iterate just once for single records
            reforg = 'ref'
            syn = Multisyn(None, {}, None)
            for sample, samplecell in zip(samples, samplecells):
                #logger.info(f"Parsing {samplecell}")
                if samplecell[i] == '-': # skip empty records
                    continue

                vals = samplecell[i].split(',')
                reforg = vals[1] # should be the same in all records, but we don't know which ones are present, so set it each time to be sure
                syn.ranges_dict[sample] = read_psf_range(vals[0], org=sample)
                
                # read cigars if present
                if len(vals) > 2:
                    if syn.cigars_dict:
                        syn.cigars_dict[sample] = cigar.cigar_from_string(vals[2])
                    else: # initialise if it hasn't been already
                        syn.cigars_dict = {sample: cigar.cigar_from_string(vals[2])}

            if reforg == 'ref':
                syn.ref = refrng
            else:
                if reforg in syn.ranges_dict:
                    syn.ref = syn.ranges_dict[reforg]
                    del syn.ranges_dict[reforg] # a ref shouldn't also be in the samples
                    # remove alignment , would just be a full match anyway
                    if syn.cigars_dict and reforg in syn.cigars_dict:
                        del syn.cigars_dict[reforg]
                else:
                    logger.error(f"Error while reading PSF: Specified reference not found in PSF!\n Line: {line}")
                    raise ValueError("Reference not found in line!")

            syns.append(syn)
    fin.close()
    return pd.DataFrame(data=list(syns)) # shouldn't require sorting

cpdef read_ancient_psf(f):
    """DEPRECATED: reads the pre-0.1 version of the PSF format. Use the new read function instead, unless working with legacy files.

    Takes a file object or path to a file in PSF format and reads it in as a DataFrame.
    """
    if isinstance(f, str):
        f = open(f, 'rt')
    syns = deque()
    orgs = f.readline().strip()[1:].split("\t")[2:] # 0 is ANN, 1 is ref
    for l in f:
        l = l.strip().split('\t')
        if l[0] == 'SYN': # line contains a multisyn region
            syn = Multisyn(read_psf_range(l[1], org="ref"), # extract reference range
                {orgs[i]:read_psf_range(cell.split(",")[0], org=orgs[i]) # extract ranges dict
                        for i, cell in enumerate(l[2:]) if cell != '.'},
                {orgs[i]:cigar.cigar_from_string(cell.split(',')[1])
                 for i, cell in enumerate(l[2:]) if cell != '.' and len(cell.split(',')) > 1} # extract cigars dict
            )
            syns.append(syn)
    f.close()
    return pd.DataFrame(data=list(syns)) # shouldn't require sorting
# END
