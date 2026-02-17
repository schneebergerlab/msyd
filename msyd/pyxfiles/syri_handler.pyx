# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3


import numpy as np
import pandas as pd

from collections import deque, defaultdict, OrderedDict

import msyd.util as util
import msyd.cigar as cigar
from msyd.multisyn import Multisyn, MultisynContainer, ChromContainer
from msyd.coords import Range
from msyd.orgs import OrgContainer

import logging
logger = util.CustomFormatter.getlogger(__name__)

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

#cpdef extract_syri_snvs(fin):
#    syri_regs = deque()
#    with open(f, 'r') as fin:
#        for line in fin:
#            l = line.strip().split()
#            if l[10] == 'SNP':
#                #TODO maybe store annotation information from fields 8-10
#                snv = SNV(Position('a', 'x', l[0], int(l[1])), Position('b', 'x', l[5], int(l[6])), l[4], l[5])
#                syri_regs.append(snv)
#
#    df = pd.DataFrame(list(syri_regs))#[[0, 1, 3, 4, 5, 6, 8, 9, 10]]
#    #TODO maybe do chromosome mapping?
#    return df

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

# given a bam file and corresponding SYNAL range df,
# Transform them into one list of Multisyn objects
cpdef match_synal(syndf, alndf, ref='a', refname="ref"):
    """
    This function takes an aligment and SYNAL dataframe and matches corresponding regions.
    It returns a dataframe containing the regions with the corresponding CIGAR string as a `Multisyn` object.
    :params: syndf: SYNAL dataframe, alndf: alignment dataframe, ref: whether the reference is the 'a' or 'b' strand in the alignment dataframe.
    :returns: a dataframe containing the SYNAL regions with corresponding CIGAR strings as `Multisyn` objects.
    """
    cdef:
        ret = MultisynContainer(None) # org info is set when returning
        syniter = syndf.iterrows()
        alniter = alndf.iterrows()
        str refchr = ref + "chr"
        str refstart = ref + "start"
        str refend = ref + "end"
        synr = next(syniter)[1]
        alnr = next(alniter)[1]
        counter = 0

    while True:
        counter += 1
        try:
            org = synr[1].org
            if synr[0].chrom == alnr[refchr] and synr[0].start == alnr[refstart] and synr[0].end == alnr[refend]:
                cg = cigar.cigar_from_string(alnr['cg'])
                rng = synr[1]
                multisyn = Multisyn(ref=synr[0], ranges_dict={org:rng}, cigars_dict={org:cg})
                #rng.end = rng.start + cg.get_len(ref=False) -1

                ## Correct mismatches between CIGAR and coordinate len

                # check on ref
                if not len(multisyn.ref) == cg.get_len():
                    logger.warning(f"CIGAR len ({cg.get_len()}) not matching coordinate len ({len(multisyn.ref)}) on ref! Adjusting end to match CIGAR len (this might be because of clipping).")
                    multisyn.ref.end = multisyn.ref.start + cg.get_len() - 1

                # check on org
                if not len(rng) == cg.get_len(ref=False):
                    # forcibly ajust end position to match cigar length, as that doesn't always seem to be the case in syri/pysam output for some reason
                    logger.warning(f"CIGAR len ({cg.get_len(ref=False)}) not matching coordinate len ({len(rng)}) on {org}! Adjusting end to match CIGAR len (this might be because of clipping).")
                    rng.end = rng.start + cg.get_len(ref=False) -1
                
                ret.append(multisyn)
                synr = next(syniter)[1]
            alnr = next(alniter)[1]
        except StopIteration:
            break

    if len(ret) <= 0.1*counter:
        logger.error("Less than 10% of syns had a matching alignment! Check that syri was run on the same alignment as was provided!")
    ret.orgs = OrgContainer.from_list([refname, org])
    return ret

cpdef handle_conflicts(syniter):
    """
    Part of the preprocessing of SYNAL regions for find_multisyn.
    Removes overlap from the first region if two overlapping regions are next to each other.
    Assumes syn to be sorted.
    Mutates syn, returns nothing.
    """
    try:
        prev = next(syniter)
    except StopIteration:
        logger.error("handle_conflicts called on empty synteny list! Most likely there is an issue with reading the input files.")

    for cur in syniter:
        #logger.debug(f"Prev: {prev}")
        #logger.debug(f"Cur: {cur}")
        if cur.ref.chrom != prev.ref.chrom: # there can be no overlap between chrs
            prev = cur
            continue

        ## check for & remove overlap on the reference
        ov = prev.ref.end - cur.ref.start +1
        if ov > 0:
            # there is overlap on ref
            logger.warning(f"Found {ov} bp overlapping synteny on {cur.ref.org} at {cur.ref.start}, trimming latter record!")
            logger.debug(f"Cur before dropping: {cur}")
            cur.drop_inplace(ov, 0) # call drop_inplace to mutate the dataframe from a reference
            logger.debug(f"Cur after dropping: {cur}")

        ## check for overlap on other orgs

        # when this is called, cur and prev should normally have the same orgs
        # will not catch overlap between non-adjacent regions!
        #assert(set(cur.ranges_dict) == set(prev.ranges_dict))
        for org in cur.ranges_dict: # should be on the same chrom
            if org not in prev.ranges_dict or cur.ranges_dict[org] is None or prev.ranges_dict[org] is None:
                continue
            assert(cur.ranges_dict[org].chrom == prev.ranges_dict[org].chrom) # prev.ranges_dict[org] is None sometimes?? O.o

            ov = prev.ranges_dict[org].end - cur.ranges_dict[org].start + 1 # indices are inclusive
            if ov > 0:
                # check if the region is fully contained, in case this ever happens
                # drop the region on this org in that case
                if cur.ranges_dict[org].end <= prev.ranges_dict[org].end:
                    logger.warning(f"On {org}, a syntenic region fully contains another! Dropping {org} from contained region.")
                    logger.debug(f"{cur.ranges_dict[org]} contained in {prev.ranges_dict[org]}!")
                    #del cur.ranges_dict[org] # this causes a crash during iteration
                    cur.ranges_dict[org] = None # set to None instead
                    if cur.cigars_dict:
                        del cur.cigars_dict[org]
                    # technically, on `org` the following syns should now be compared
                    # however that would require storing the last region for every org separately
                    # for a case that shouldn't ever occur
                    # => just delete it and skip comparisons
                    continue

                # there is overlap on org
                logger.warning(f"Found {ov} bp overlapping synteny on {org} at {cur.ranges_dict[org].start}, trimming latter record!")
                logger.debug(f"Overlapping on {org}: {prev}, {cur}")
                logger.debug(f"Cur before dropping: {cur}")
                cur.drop_on_org_inplace(ov, 0, org)
                logger.debug(f"Cur after dropping: {cur}")

        prev = cur
# END
