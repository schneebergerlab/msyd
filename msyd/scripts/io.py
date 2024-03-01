from collections import deque, OrderedDict
import logging
import numpy as np
import pandas as pd
from msyd.coords import Range
VARS = ['SYN', 'SYNAL', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']
def readsyriout(f):
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

def extract_syri_regions(rawsyriout, ref='a', anns=['SYN'], reforg='ref', qryorg='qry'):
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
    merged = pd.concat([rawsyriout.loc[rawsyriout['type'] == ann if 'type' in rawsyriout.columns else rawsyriout['vartype'] == ann] for ann in anns]) # different syri versions seem to use different names for the type
    # if implementing filtering later, filter here

    for row in merged.iterrows():
        row = row[1]
        # removed util.chrom_to_int, was causing problems
        buf.append([Range(reforg, row[refchr], refhaplo, row[refstart], row[refend]),
            Range(qryorg, row[qrychr], qryhaplo, row[qrystart], row[qryend])
            ])

    return pd.DataFrame(data=list(buf), columns=[reforg, qryorg])


def extract_syri_regions_from_file(fin, ref='a', anns=['SYN'], reforg='ref', qryorg='qry'):
    raw, chr_mapping = readsyriout(fin) #TODO? handle chr_mapping
    return extract_syri_regions(raw, ref=ref, anns=anns, reforg=reforg, qryorg=qryorg)


def extract_syri_regions_to_list_from_files(fins, qrynames, cores=1, **kwargs):
    """
    `extract_syri_regions`, but for processing a list of inputs
    """
    if len(fins) != len(qrynames):
        logger.error(f"Infiles and qrynames lists lengths not matching. Offending lists: {fins} and {qrynames}")
    partial = lambda x, qryname: extract_syri_regions_from_file(x, qryorg=qryname, **kwargs)

    if cores == 1:
        syns = [partial(fin, qryname) for fin, qryname in zip(fins, qrynames)]
    else:
        # `partial` requires two parameters, only 1 is given here. would crash ?
        with Pool(cores) as pool:
            syns = pool.map(partial, fins)

    return syns
    #return [extract_syri_regions(fin, **kwargs,\
    #        #reforg=fin.split('/')[-1].split('_')[0],\
    #        qryorg=fin.split('/')[-1].split('_')[-1].split('syri')[0])\
    #        for fin in fins]
# END