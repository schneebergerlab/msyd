#!/bin/bash

for seq in seqs/*.fna.gz
do
	bs=$(basename $seq)
	bsub -q multicore20 -n5 -R"span[hosts=1] rusage[mem=2000]" -M10000 \
	-oo job_syri_${bs}.log -eo job_syri_${bs}.err \
	"syri --nc 5 -F P --cigar --dir syri --prefix $bs -c alns/$bs.paf -r ref.fna.gz -q $seq --lf $bs.syri.log"
done
