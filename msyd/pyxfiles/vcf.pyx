# -*- coding: utf-8 -*-
# distutils: language = c++
# cython: language_level = 3


import os
import re
from copy import copy

import pysam

#from msyd.vars import SNV
import msyd.util as util
import msyd.io as io

logger = util.CustomFormatter.getlogger(__name__)

HEADER="""##INFO=<ID=END,Number=1,Type=Integer,Description="End position on reference genome">
##ALT<ID=CORESYN,Description="Coresyntenic region (syntenic between any two samples)">
##ALT<ID=MERASYN,Description="Merasyntenic region (syntenic between any two samples for a strict subset of the samples)">
##INFO=<ID=PID,Number=1,Type=Integer,Description="Numerical part of the ID of the parent PANSYN region. If the PID of a region is 10, it's parent's ID will be CROSSSYN10 or CORESYN10 (and there will be only one of either).">
##FORMAT=<ID=CHR,Number=1,Type=String,Description="Chromosome in this sample">
##FORMAT=<ID=START,Number=1,Type=Integer,Description="Start position in this sample">
##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position  in this sample">
##FORMAT=<ID=CG,Number=1,Type=String,Description="CIGAR String containing the alignment to the reference">
##FORMAT=<ID=AI,Number=1,Type=Integer,Description="Alignment Identity of the alignment of this region to the reference">
##FORMAT=<ID=SYN,Number=1,Type=Integer,Description="1 if this region is syntenic to reference, else 0">"""
##FORMAT=<ID=HAP,Number=1,Type=Character,Description="Unique haplotype identifier">"""

cpdef filter_vcfs(syns, vcfs: List[Union[str, os.PathLike]], ref: Union[str, os.PathLike], add_syn_anns=False, no_complex=False, impute_ref=False):
    tmpfiles = [util.gettmpfile() for _ in vcfs]

    for i in range(len(vcfs)):
        logger.info(f"Filtering {vcfs[i]}")
        syri_vcf = not re.fullmatch(r".*syri\.vcf$", vcfs[i]) == None
        extract_syntenic_from_vcf(syns, vcfs[i], tmpfiles[i], ref=ref, add_syn_anns=add_syn_anns, no_complex=no_complex, coords_in_info=syri_vcf, impute_ref=impute_ref)

    return tmpfiles

cpdef void extract_syntenic_from_vcf(syns, inpath:Union[str, os.PathLike], outpath: Union[str, os.PathLike], force_index=True, synorg='ref', ref=None, keep_anns=True, add_syn_anns=True, add_cigar=False, add_identity=True, no_complex=False, coords_in_info=False, impute_ref=False):
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
        int syncounter = 1

    if not set(vcfin.header.samples).issubset(orgs):
        removing = set(vcfin.header.samples).difference(orgs)
        logger.warning(f"Input VCF contains organisms not in PSF file! Double-Check names used in .tsv. Removing samples {removing} from VCF")
        vcfin.subset_samples(orgs)

    # read reference if it hasn't been read already
    if ref and type(ref) != dict:
        logger.info("Reading in Reference Fasta")
        ref = io.readfasta(ref)
    elif not ref:
        logger.warning("No Reference specified, not saving Ref Sequence in VCF!")

    # add header required for storing PANSYN annotations
    if add_syn_anns or coords_in_info:
        for line in HEADER.splitlines():
            vcfout.header.add_line(line)

    orgsvcf = list(vcfin.header.samples) # select only contained organisms

    if coords_in_info and len(orgsvcf) != 1:
        logger.error("reading coords from INFO only supported for VCFs with exactly one sample! Check if your SyRI installation is up to date!")
        raise ValueError("reading coords from INFO only supported for VCFs with exactly one sample!")

    # force indexing to allow for calling fetch later.
    #TODO try using until_eof=True as mentioned in the pysam FAQ
    if force_index and not vcfin.index_filename:
        vcfin.close()
        pysam.tabix_index(inpath, force=True, preset='vcf', keep_original=True)
        inpath += ".gz" # no way to turn off automatic compression, apparently
        vcfin = pysam.VariantFile(inpath)

    # add multisyn regions and contained records
    for syn in syns.iterrows():
        syn = syn[1][0]
        rng = syn.ref
        if synorg != 'ref':
            # untested as of now
            if synorg in syn.ranges_dict:
                rng = syn.ranges_dict[synorg]
                syn = copy(syn)
                syn.ref = rng
                del syn.ranges_dict[synorg]
            else: # ignore regions not present in this org
                continue
        
        # add the multisyn region, if specified
        if add_syn_anns:
            try:
                add_syn_ann(syn, vcfout, ref=ref, no=syncounter, add_cigar=add_cigar, add_identity=add_identity)
            except ValueError:
                logger.error(f"Error adding multisyn annotation for region {syn} to VCF. Check if the chromosome names match!")
            syncounter +=1

        # debugging
        #pos = 0
        #oldrec = None

        # write the small variants in the multisyn region
        for rec in vcfin.fetch(rng.chr, rng.start, rng.end + 1): # pysam is half-inclusive
            # double check if the chr has been added, was throwing errors for some reason...
            if rec.chrom not in header_chrs:
                #logger.info(f"extract_from_syntenic Adding {rec.chrom} to header")
                header_chrs.add(rec.chrom)
                if ref:
                    # add length if it is known from the reference
                    vcfout.header.add_line("##contig=<ID={},length={}>".format(rec.chrom, len(ref[rec.chrom])))
                else:
                    vcfout.header.add_line("##contig=<ID={}>".format(rec.chrom))

            # do not include vars that start before this region
            # otherwise the sorting will be violated
            if rec.pos < rng.start:
                continue

            # check if the region is complex by looking for symbolic alleles
            if no_complex:
                #if any(map(lambda x: re.fullmatch(r'N|[ACGT]*', x) == None, rec.alleles)):
                # does this need checking for None alleles? not sure...
                if any([re.fullmatch(r'N|[ACGT]*', allele) == None for allele in rec.alleles]):
                    logger.debug(f"Complex variant: skipping {rec}")
                    continue # skip this variant

            # debugging
            #if pos > rec.pos and rec.chrom == oldrec.chrom:
            #    logger.error(f"Unsorted! {oldrec} before {rec}! ({pos} > {rec.pos}")
            #    pos = rec.pos
            #    oldrec = rec
            #logger.debug(f"syn: {rng.start}-{rng.end}, pos {rec.pos}")

            new_rec = vcfout.new_record()
            new_rec.pos = rec.pos
            new_rec.chrom = rec.chrom
            new_rec.id = rec.id
            new_rec.alleles = rec.alleles
            #new_rec.ref = rec.ref # unnecessary, as this is covered by alleles

            # discard old INFO information if reading it in as coords
            if not coords_in_info:
                for key in rec.info:
                    new_rec.info[key] = rec.info[key]
            # add Parent information
            if add_syn_anns:
                new_rec.info['PID'] = syncounter
            
            ## Handle the sample fields
            for sample in rec.samples:
                # see also
                # https://github.com/pysam-developers/pysam/blob/43c10664834cf3914000a744e6057826c6a6fa65/pysam/libcbcf.pyx#L3443
                # this stuff is all undocumented (see pysam#407)

                # transfer over existing annotations if instructed
                if keep_anns:
                    new_rec.samples[sample].update(rec.samples[sample])

                # impute ref if specified and syntenic to ref
                if impute_ref\
                        and sample in syn.ranges_dict \
                        and rec.samples[sample].alleles[0] is None:
                    #logger.info(f"imputing {sample}")
                    # no genotype information present
                    # annotate as reference genotpye if impute is enabled
                    if syn.ref.org == 'ref':
                        new_rec.samples[sample].allele_indices = 0
                    else: # should catch non-ref. merasynteny
                        new_rec.samples[sample].update(rec.samples[syn.ref.org])
                        #new_rec.samples[sample].allele_indices = rec.samples[syn.ref.org].allele_indices

            # read in coords from INFO column, add to single sample
            # TODO get this to work, also re-look at the if below, seeems not right (rec shouldn't be writeable)
            if coords_in_info:
                sample = orgsvcf[0] # there can only be one sample
                for info, ft in [('ChrB', 'CHR'), ('StartB', 'START'), ('EndB', 'END')]:
                    if info in rec.info:
                        new_rec.samples[sample][ft] = rec.info[info]

            vcfout.write(new_rec)
    #vcfout.close()
    #vcfin.close()

cpdef void reduce_vcfs(vcfs: List[Union[str, os.PathLike]], opath: Union[str, os.PathLike], add_syn_anns=True):
    # quick and dirty reduction function, TODO write proper one when integrating with PSF variation merging

    if len(vcfs) < 1:
        logger.error("reduce_vcfs called with empty vcfs!")
        return
    elif len(vcfs) == 1:
        logger.warning(f"reduce_vcfs called with only one vcf: {vcfs}")
        return
    elif len(vcfs) == 2:
        merge_vcfs(vcfs[0], vcfs[1], opath)
        return

    tmpfiles = [util.gettmpfile() for _ in range(2, len(vcfs))] # the first two and last vcfs don't need to be stored as tempfiles
    merge_vcfs(vcfs[0], vcfs[1], tmpfiles[0])
    for i in range(1, len(vcfs)-2):
        merge_vcfs(tmpfiles[i-1], vcfs[i+1], tmpfiles[i])
    # incorporate the last vcf, save directly to output
    merge_vcfs(vcfs[-1], tmpfiles[-1], opath)

#cpdef void impute_syn_in_vcf(syns, inpath: Union[str, os.PathLike], outpath: Union[str, os.PathLike], force_index=True, synorg='ref', ref=None, keep_anns=True, add_syn_anns=True, add_cigar=False, add_identity=True, no_complex=False, coords_in_info=False, impute_ref=False):



cpdef add_syn_anns_to_vcf(syns, vcfin: Union[str, os.PathLike], vcfout: Union[str, os.PathLike], ref=None):
    """Takes a VCF file, overwrites it adding annotations for core/cross-syn region. Other records are preserved as-is."""
    cdef:
        # read in old file, initialise for overwriting
        oldvcf = pysam.VariantFile(vcfin, 'r')
        newvcf = pysam.VariantFile(vcfout, 'w')
        int syncounter = 1

    # copy header, deduplicate along the way
    headerset = set()
    for line in str(oldvcf.header).splitlines()[:-1]:
        #logger.info(line)
        if not line in headerset:
            newvcf.header.add_line(line)
            headerset.add(line)

    # extend header if necessary
    for line in HEADER.splitlines():
        if not line in headerset:
            newvcf.header.add_line(line)
            headerset.add(line)

    for sample in oldvcf.header.samples:
        newvcf.header.add_sample(sample)

    syniter = syns.iterrows()
    syn = next(syniter)[1][0]

    for oldrec in oldvcf:
        if oldrec.chrom < syn.ref.chr:
            if oldrec.chrom not in set(newvcf.header.contigs): # check if chr needs adding
                #logger.info(f"add_syn_anns_to_vcf Adding {syn.ref.chr} to header")
                newvcf.header.add_line("##contig=<ID={}>".format(oldrec.chrom))
            copy_record(oldrec, newvcf, pid=syncounter-1)
        elif oldrec.start < syn.ref.start:
            copy_record(oldrec, newvcf, pid=syncounter-1)
        else:
            # the record is not before the syn
            add_syn_ann(syn, newvcf, ref=ref, no=syncounter)
            copy_record(oldrec, newvcf, pid=syncounter)
            syncounter += 1
            try:
                syn = next(syniter)[1][0]
            except StopIteration:
                syn.ref.start = 99999999999 # to disable ever entering the else case again


cdef add_syn_ann(syn, ovcf, ref=None, no=None, add_cigar=False, add_identity=True):
    rng = syn.ref
    rec = ovcf.new_record()
    rec.start = rng.start
    rec.pos = rec.start
    rec.stop = rng.end
    chrom = rng.chr

    if chrom not in set(ovcf.header.contigs):
        #logger.info(f"add_syn_ann Adding {chrom} to header")
        if ref:
            # add length if it is known from the reference
            ovcf.header.add_line("##contig=<ID={},length={}>".format(chrom, len(ref[chrom])))
        else:
            ovcf.header.add_line("##contig=<ID={}>".format(chrom))

    rec.chrom = chrom
    if set(ovcf.header.samples).issubset(syn.get_orgs()): # if the region is coresyn within the VCF
        if ref:
            rec.alleles = [ref[rec.chrom][rec.start], "<CORESYN>"]
        else:
            rec.alleles = ["<SYN>", "<CORESYN>"]
        rec.id = "CORESYN{}".format(no)
    else:
        if ref:
            rec.alleles = [ref[rec.chrom][rec.start], "<CROSSSYN>"]
        else:
            rec.alleles = ["<SYN>", "<CROSSSYN>"]
        rec.id = "CROSSSYN{}".format(no)

    #rec.info['NS'] = syn.get_degree() # update NS column, include not only orgs in sample now

    # write the multisyn annotation
    for org in ovcf.header.samples:
        if org in syn.get_orgs():
            rng = syn.ranges_dict[org]
            rec.samples[org].update({'SYN':1, 'CHR':rng.chr, 'START': rng.start, 'END': rng.end})
            if syn.cigars_dict:
                cg = syn.cigars_dict[org]
                if add_cigar:
                    rec.samples[org].update({'CG': cg.to_string()})
                if add_identity:
                    rec.samples[org].update({'AI': int(cg.get_identity()*100)})
        else:
            rec.samples[org].update({'SYN': 0})
    # write to file
    ovcf.write(rec)


cdef str merge_vcfs(lf: Union[str, os.PathLike], rf:Union[str, os.PathLike], of:Union[str, os.PathLike]):
    logger.info(f"Merging {lf} and {rf} into {of}")
    # TODO reimplement this with common framework with merge psfs
    # do all this in memory to be faster
    lvcf = pysam.VariantFile(lf, 'r')
    rvcf = pysam.VariantFile(rf, 'r')
    ovcf = pysam.VariantFile(of, 'w')

    # Prepare the header
    if str(lvcf.header) != str(ovcf.header):
        logger.info(f"Headers not matching in {lf} and {rf}! Combining.")

    # merge the headers, deduplicate along the way
    headerset = set()
    for line in str(lvcf.header).splitlines()[1:-1]:
        if not line in headerset:
            ovcf.header.add_line(line)
    for line in str(rvcf.header).splitlines()[1:-1]:
        if not line in headerset:
            ovcf.header.add_line(line)

    # Panic on two empty vcfs
    if len(rvcf.header.samples) == 0 or len(lvcf.header.samples) == 0:
        logger.error("Merging VCFs with no samples is not supported, exiting!")
        #return

    logger.info(f"Found samples: {list(lvcf.header.samples)}, {list(rvcf.header.samples)}")
    for sample in rvcf.header.samples:
        if sample in ovcf.header.samples:
            logger.warning(f"Duplicate sample '{sample}' encountered in {rf}. Appending filename to sample name.")
            sample += '_'
            sample += rf
        ovcf.header.add_sample(sample)
    for sample in lvcf.header.samples:
        if sample in ovcf.header.samples:
            logger.warning(f"Duplicate sample '{sample}' encountered in {lf}. Appending filename to sample name.")
            sample += '_'
            sample += lf
        ovcf.header.add_sample(sample)
    #print(ovcf.header)


    # iterate through the VCFs, merging individual records
    # keeps only records that can be merged
    lann = None
    rann = None
    try:
        lann = next(lvcf)
        rann = next(rvcf)
    except StopIteration:
        logger.error("Empty VCF encountered. Outputting empty VCF!")
        return of

    try:
        while True:
            # skip until we are at the same position
            if lann.chrom != rann.chrom:
                 # check if chr matches, otherwise skip till end
                if lann.chrom < rann.chrom:
                    if lann.chrom not in set(ovcf.header.contigs):
                        #logger.info(f"merge_vcfs Adding {lann.chrom} to header")
                        ovcf.header.add_line("##contig=<ID={}>".format(lann.chrom))
                    # skip till end
                    while lann.chrom < rann.chrom:
                        copy_record(lann, ovcf)
                        lann = next(lvcf)
                    continue
                else:
                    if rann.chrom not in set(ovcf.header.contigs):
                        logger.info(f"merge_vcfs Adding {rann.chrom} to header")
                        ovcf.header.add_line("##contig=<ID={}>".format(rann.chrom))
                    # skip till end
                    while rann.chrom < lann.chrom:
                        copy_record(rann, ovcf)
                        rann = next(rvcf)
                    continue

            elif lann.pos < rann.pos:
                copy_record(lann, ovcf)
                lann = next(lvcf)
                continue
            elif rann.pos < lann.pos or rann.chrom < lann.chrom:
                copy_record(rann, ovcf)
                rann = next(rvcf)
                continue
            # lann.pos and rann.pos must be equal
            elif lann.alleles[0] != rann.alleles[0]:
                # there are multiple annotations on this position, and they aren't sorted
                # (if they are sorted the code above already works)
                pos = lann.pos
                chrom = lann.chrom

                # extract all records at this position in rann, store in dictionary
                mapping = {rann.alleles[0]: rann} if rann.alleles[0] != 'N' else dict()
                rann = next(rvcf)
                while rann.pos == pos and rann.chrom == chrom:
                    if rann.alleles[0] == 'N': # to handle NOTAL and HDR records
                        pass
                    #elif rann.alleles[0] in mapping:
                    #    logger.error(f"Two variants with same reference at same position in same file({rf}): {rann}, {mapping[rann.alleles[0]]}.")
                    #    #merge_vcf_records(rann, mapping[rann.alleles[0]], ovcf)
                    else:
                        mapping[rann.alleles[0]] = rann
                    rann = next(rvcf)
                # rann now contains the first record after this position

                # match with records in lvcf one by one
                while lann.pos == pos and lann.chrom == chrom:
                    if lann.alleles[0] != 'N' and lann.alleles[0] in mapping:
                        merge_vcf_records(lann, mapping[lann.alleles[0]], ovcf)
                        del mapping[lann.alleles[0]]
                    else:
                        copy_record(lann, ovcf)
                    lann = next(lvcf)
                # lann now contains the first record after this position

                # add records in rvcf that do not match any in lvcf
                for record in mapping.values():
                    copy_record(record, ovcf)

                # all variants up to this position have been added, continue as normal
                continue
            else:
                # they are at the same position, and the ref gt matches
                merge_vcf_records(lann, rann, ovcf)
                # discard these two records, look at the next
                lann = next(lvcf)
                rann = next(rvcf)
                continue

    except StopIteration:
        pass

    return of # to enable reduction operation

cdef copy_record(rec: VariantRecord, ovcf:VariantFile, int pid=0):
    """
    Utility function to copy a record to another VCF, because pysam needs some conversions done.
    """
    if rec.chrom not in set(ovcf.header.contigs):
        logger.info(f"copy_record Adding {rec.chrom} to header")
        ovcf.header.add_line("##contig=<ID={}>".format(rec.chrom))
    new_rec = ovcf.new_record()
    new_rec.pos = rec.pos
    new_rec.chrom = rec.chrom
    new_rec.id = rec.id
    new_rec.alleles = rec.alleles
    for key in rec.info:
        new_rec.info[key] = rec.info[key]
    for sample in rec.samples:
        new_rec.samples[sample].update(rec.samples[sample])
    if pid != 0: # set parent ID if necessary
        new_rec.info['PID'] = pid
    ovcf.write(new_rec)

cdef merge_vcf_records(lrec: VariantRecord, rrec:VariantRecord, ovcf:VariantFile):
    """
    Merge two vcf records from different files, append to ovcf.
    """
    rec = ovcf.new_record()
    rec.pos = lrec.pos # shoud be equal anyway

    chrom = lrec.chrom
    # this should not be necessary, but for some reason the chrs do not seem to be added by merging the header?
    if chrom not in set(ovcf.header.contigs):
        logger.info(f"merge_vcf_records Adding {chrom} to header")
        ovcf.header.add_line("##contig=<ID={}>".format(chrom))

    rec.chrom = chrom

    lref = lrec.alleles[0]
    rref = rrec.alleles[0]

    rec.stop = lrec.stop

    if lref != rref:
        logger.error("Trying to join records with different references!")

    # construct joined gt -> index map
    gtmap = {gt:ind+1 for ind, gt in enumerate(set(lrec.alleles[1:] + rrec.alleles[1:]))}

    alleles = [rref] + list(gtmap.keys())
    gtmap[rref] = None # add the reference to gtmap

    # <NOTAL> annotations have only one allele in SyRI VCF files
    # pysam throws an error when storing variants with only one allele,
    # but can read them just fine
    if len(alleles) == 1:
        alleles.append(' ') # try to trick pysam
    rec.alleles = alleles

    if lrec.id != rrec.id:
        logger.warning(f"id not matching in {lrec.id} and {rrec.id}! Choosing {lrec.id}")
    rec.id = lrec.id


    # rec.samples.update() throws internal pysam errors, circumvent it by directly calling update for each sample
    for samples in [lrec.samples, rrec.samples]:
        for sample in samples:
            rec.samples[sample].update(samples[sample])

    # handle GT column separately, incorporating the gtmap constructed earlier
    #for samples in [lrec.samples, rrec.samples]:
    for (samples, alleles) in [(lrec.samples, lrec.alleles), (rrec.samples, rrec.alleles)]:
        for sample in samples:
            if not 'GT' in rec.samples[sample] or alleles == None: # nothing needs updating
                continue
            # apparently pysam treats the genotype specially without documenting that behaviour...
            gt = rec.samples[sample]['GT']
            mapper = lambda x: gtmap[alleles[x]] if x is not None else gtmap[alleles[0]]
            if not gt:
                logger.warning(f"Invalid GT found: {gt} for {sample} in {rec.id}")
                continue
            elif len(gt) == 2:
                rec.samples[sample]['GT'] = (mapper(gt[0]), mapper(gt[1]))
            elif len(gt) == 1:
                # there is an unphased GT
                rec.samples[sample]['GT'] = mapper(gt[0])
            else:
                # there is an invalid GT
                logger.warning(f"Invalid GT found: {gt} for {sample} in {rec.id}")

            #logger.info(f"{gt}, {gtmap}, {alleles}, {rec.samples[sample]['GT']}")

    #if list(lrec.format) != list(rrec.format):
    #    # temporary prints necessary because pysam is annoying
    #    logger.info(f"format not matching: {list(lrec.format)} and {list(rrec.format)}!")

    # pysam does not allow setting the info field all at once, do it iteratively:
    for key in lrec.info:
        rec.info[key] = lrec.info[key]

    for key in rrec.info:
        if key in lrec.info and lrec.info[key] != rrec.info[key]:
            logger.warning(f"Conflicting info stored for {key} in {rec.id}: {lrec.info[key]} != {rrec.info[key]}! Choosing {lrec.info[key]}")
            #continue
        else:
            rec.info[key] = rrec.info[key]

    ovcf.write(rec)

