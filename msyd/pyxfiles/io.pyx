# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3

import numpy as np
import pandas as pd
from scipy.stats import *
import pysam

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

from msyd.coords import Range
from msyd.multisyn import Multisyn
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
                syri_regs.append(snv)

    df = pd.DataFrame(list(syri_regs))#[[0, 1, 3, 4, 5, 6, 8, 9, 10]]
    #TODO maybe do chromosome mapping?
    return df

# cython-lint flags the default arg list as dangerous
# but in this case it's fine since its static
cpdef extract_syri_regions_from_file(fin, ref='a', anns=['SYN'], reforg='ref', qryorg='qry'): # no-cython-lint
    raw, _chr_mapping = readsyriout(fin) #TODO? handle chr_mapping
    return extract_syri_regions(raw, ref=ref, anns=anns, reforg=reforg, qryorg=qryorg)


# cython-lint flags the default arg list as dangerous
# but in this case it's fine since its static
cpdef extract_syri_regions(rawsyriout, ref='a', anns=['SYN'], reforg='ref', qryorg='qry'): # no-cython-lint
    """
    Given a syri output file, extract all regions matching a given annotation.
    Returns the output as a dict containing one Dataframe per chromosome.
    """
    # columns to look for as start/end positions
    refchr = ref + "chr"
    refstart = ref + "start"
    refend = ref + "end"

    qry = 'b' if ref == 'a' else 'a' # these seem to be the only two values in syri output
    qrychr = qry + "chr"
    qrystart = qry + "start"
    qryend = qry + "end"


    merged = pd.concat([rawsyriout.loc[rawsyriout['type'] == ann if 'type' in rawsyriout.columns else rawsyriout['vartype'] == ann] for ann in anns]) # different syri versions seem to use different names for the type
    if merged.empty:
        logger.error(f"No annotation of type in {anns} found!")

    out = dict()
    buf = deque()
    chrom = merged.iloc[0].at[refchr] #merged.at[1, refchr] # throws an error if the first index is not 1
    for _, row in merged.iterrows():
        # write buffer to out if necessary
        if row[refchr] != chrom:
            out[chrom] = pd.DataFrame(data=list(buf), columns=[reforg, qryorg])
            chrom = row[refchr]
            buf = deque()
        # append current line to buffer
        buf.append([Range(reforg, row[refchr], row[refstart], row[refend]),
            Range(qryorg, row[qrychr],  row[qrystart], row[qryend])
            ])

    # add last chr
    out[chrom] = pd.DataFrame(data=list(buf), columns=[reforg, qryorg])

    return out

def extract_from_filelist(fins, qrynames, cores=1, **kwargs):
    """
    `extract_syri_regions`, but for processing a list of inputs.
    Will return a 
    """
    if len(fins) != len(qrynames):
        logger.error(f"Infiles and qrynames lists lengths not matching. Offending lists: {fins} and {qrynames}")

    out = defaultdict(list)
    # optionally parallelize i/o like this?
#    with Pool(cores) as pool:
#        annoying_workaround = partial(extract_syri_regions_from_file, **kwargs)
#        for chrom, syndf in pool.map(annoying_workaround, zip(fins, qrynames)):
    for fin, qryname in zip(fins, qrynames):
        for chrom, syndf in extract_syri_regions_from_file(fin, qryorg=qryname, **kwargs).items():
            out[chrom].append(syndf)

    return out

cpdef void save_to_vcf(syns: Union[str, os.PathLike], outf: Union[str, os.PathLike], ref=None, cores=1, add_cigar=False, add_identity=True):
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

    if ref and type(ref) != dict:
        logger.info("Reading in Reference Fasta")
        ref = readfasta(ref)
    elif not ref:
        logger.warning("No Reference specified, not saving Ref Sequence in VCF!")

    #out.header.add_samples(util.get_orgs_from_df(syns)) # according to the documentation, this works, but the function doesn't seem to exist...
    for org in orgs:
        out.header.add_sample(org)

    # add each multisyn object
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
            rec.id = "CORESYN{}".format(corecounter)
            corecounter += 1
        else:
            if ref:
                rec.alleles = [ref[rec.chrom][rec.start], "<CROSSSYN>"]
            else:
                rec.alleles = ["<SYN>", "<CROSSSYN>"]
            rec.id = "CROSSSYN{}".format(crosscounter)
            crosscounter += 1

        #rec.info['NS'] = syn.get_degree() # update NS column, include not only orgs in sample now

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

cpdef save_to_psf(dfmap, buf, save_cigars=True, force_ref_pos=False):
    """
    Takes a map of chrom IDs to DFs containing multisyns and writes them to buf.
    Preserves the sorting of the DFs, sorts chroms lexicallicaly.
    Calls to `save_df_to_psf`.
    """
    if len(dfmap) == 0:
        raise ValueError("Empty dfmap provided!")

    # write header; assumes the first chrom contains all orgs at least once
    buf.write("#CHR\tSTART\tEND\tANN\tREP\tRCHR\tRSTART\tREND\t")
    buf.write("\t".join(util.get_orgs_from_df(list(dfmap.values())[0])))
    buf.write("\n")

    # write contents
    #TODO parallelize?
    for chrom in sorted(dfmap):
        save_df_to_psf(dfmap[chrom], buf, emit_header=False, save_cigars=save_cigars, force_ref_pos=force_ref_pos)

cpdef save_df_to_psf(df, buf, save_cigars=True, emit_header=True, force_ref_pos=False):
    """Takes a df containing `Multisyn` objects and writes them in population synteny file format to `buf`.
    Can be used to print directly to a file, or to print or further process the output.
    """
    # output organisms in lexicalic ordering
    orgs = sorted(util.get_orgs_from_df(df))
    cdef:
        int n = len(orgs) + 1 # to account for ref
        int corecounter = 1
        int meracounter = 1
        int privcounter = 1
        int coreend = 0
        str corechr = ''

    if emit_header:
        buf.write("#CHR\tSTART\tEND\tANN\tREP\tRCHR\tRSTART\tREND\t")
        buf.write("\t".join(orgs))
        buf.write("\n")

    syniter = df.iterrows()
    # TODO: assert that the columns and columns are in same order as the input file (genomes.csv)
    while True:
        mesyns = []
        refmesyns = []
        privs = [] # ref private is handled during writing
        syn = None
        try:
            syn = next(syniter)[1][0]
            # get all mesyns, separate by those having a position on reference and those that don't
            while syn.get_degree() < n:
                if syn.ref.org == "ref":
                    refmesyns.append(syn)
                else:
                    if syn.get_degree() > 1:
                        mesyns.append(syn)
                    else:
                        privs.append(syn)
                syn = next(syniter)[1][0]
        except StopIteration: # try/catch block internal, so things still get written after we run out of multisyn regions
            pass

        # first, write non-ref-position merasynteny
        # write to the first position it can be
        # maybe this should be annotated for the entire range it can be instead (coreend+1:syn.start-1)
        for mesyn in mesyns:
            if force_ref_pos:
                buf.write('\t'.join([corechr, str(coreend+1), str(coreend+1), f"MERASYN{meracounter}", mesyn.ref.org, mesyn.ref.chr, str(mesyn.ref.start), str(mesyn.ref.end), '']))
            else:
                buf.write('\t'.join(['.', '.', '.', f"MERASYN{meracounter}", mesyn.ref.org, mesyn.ref.chr, str(mesyn.ref.start), str(mesyn.ref.end), '']))
            write_multisyn(mesyn, buf, orgs, save_cigars=save_cigars)
            meracounter += 1

        for priv in privs:
            if force_ref_pos:
                buf.write('\t'.join([corechr, str(coreend+1), str(coreend+1), f"PRIVATE{privcounter}", priv.ref.org, priv.ref.chr, str(priv.ref.start), str(priv.ref.end), '']))
            else:
                buf.write('\t'.join(['.', '.', '.', f"PRIVATE{privcounter}", priv.ref.org, priv.ref.chr, str(priv.ref.start), str(priv.ref.end), '']))

            write_multisyns([priv], buf, orgs, save_cigars=save_cigars)
            privcounter += 1

        # write mesyn regions that have a position on reference at their appropriate position
        for refmesyn in refmesyns:
            ref = refmesyn.ref
            buf.write('\t'.join([ref.chr, str(ref.start), str(ref.end), f"MERASYN{meracounter}" if refmesyn.get_degree() > 1 else f"PRIVATE{privcounter}", ref.org, '.', '.', '.', '']))
            write_multisyn(refmesyn, buf, orgs, save_cigars=save_cigars)
            if refmesyn.get_degree() > 1:
                meracounter += 1
            else:
                privcounter += 1

        # write coresyn region
        if syn:
            ref = syn.ref
            coreend = ref.end
            corechr = ref.chr
            buf.write('\t'.join([ref.chr, str(ref.start), str(ref.end), f"CORESYN{corecounter}", ref.org, '.', '.', '.', '']))
            write_multisyn(syn, buf, orgs, save_cigars=save_cigars)
            corecounter += 1
        else:
            break

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
         if org in multisyn.ranges_dict else
         (multisyn.ref.to_psf() if multisyn.ref.org == org else '.') # if there is no synteny, put a .
         for org in orgs])
              )
    buf.write("\n")

cpdef read_psf(fin):
    """
    Takes a file object or path to a file in PSF format and reads it in as a DataFrame of Multisynteny objects.
    Supports the new version of PSF format; for legacy files, use the deprecated version of this function.
    """
    syns = deque()
    if isinstance(fin, str):
        fin = open(fin, 'rt')

    #CHR  START  END  ANN  REF  CHR  START  END  G1  G2  G3...
    line = fin.readline().strip().split()
    samples = line[8:] # store sample name order as in file
    # should be lexicalic, but who knows
    for line in fin:
        line = line.strip().split()
        if line == []: continue

        reforg = line[4]

        refrng = Range('ref', line[0], int(line[1]), int(line[2])) if reforg == 'ref'\
            else Range(reforg, line[5], int(line[6]), int(line[7]))

        syn = Multisyn(refrng, {}, None)

        for org, entry in zip(samples, line[8:]):
            if entry == '.' or org == reforg: # skip empty records and ref
                continue

            vals = entry.split(',')
            syn.ranges_dict[org] = read_psf_range(org, vals[0])
            
            # read cigars if present
            if len(vals) > 1:
                if syn.cigars_dict:
                    syn.cigars_dict[org] = cigar.cigar_from_string(vals[1])
                else: # initialise if it hasn't been already
                    syn.cigars_dict = {org: cigar.cigar_from_string(vals[1])}
        # add read in syn to output
        syns.append(syn)
    # clean up, return
    fin.close()
    return pd.DataFrame(data=list(syns)) # shouldn't require sorting

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
                syn.ranges_dict[sample] = read_psf_range(sample, vals[0])
                
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
            syn = Multisyn(read_psf_range("ref", l[1]), # extract reference range
                {orgs[i]:read_psf_range(orgs[i], cell.split(",")[0]) # extract ranges dict
                        for i, cell in enumerate(l[2:]) if cell != '.'},
                {orgs[i]:cigar.cigar_from_string(cell.split(',')[1])
                 for i, cell in enumerate(l[2:]) if cell != '.' and len(cell.split(',')) > 1} # extract cigars dict
            )
            syns.append(syn)
    f.close()
    return pd.DataFrame(data=list(syns)) # shouldn't require sorting
# END
