#!/bin/bash

for seq in scaffolds/*.fna.gz
do
	bs=$(basename $seq)
	bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=6000]" -M10000 \
	-oo ${tmpdir}syri_$ref$org.log -eo ${tmpdir}syri_$ref$org.err \
	"syri --nc 20 -F P --cigar --prefix syri/$bs -c alns/$bs.paf -r grch38.v14_genomic.fna.gz -q $seq --lf $bs.syri.log"
done
