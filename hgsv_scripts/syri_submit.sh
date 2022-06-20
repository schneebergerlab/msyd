#!/bin/bash

for seq in scaffolds_filtered/*.fna.bgz
do
	bs=$(basename $seq)
	bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=6000]" -M10000 \
	-oo job_syri_${bs}.log -eo job_syri_${bs}.err \
	"syri --nc 20 -F P --cigar --dir syri --prefix $bs -c alns_filtered/$bs.paf -r ref-nomtsex.fna.bgz -q $seq --lf $bs.syri.log"
done
