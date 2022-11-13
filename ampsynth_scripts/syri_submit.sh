#!/bin/bash

for seq in seqs/*.fna
do
	bs=$(basename -s .fna $seq)
	# do it locally for now, cluster is busy
	#bsub -q multicore20 -n5 -R"span[hosts=1] rusage[mem=2000]" -M10000 \
	#-oo job_syri_${bs}.log -eo job_syri_${bs}.err \
	#"syri --nc 5 -F P --cigar --dir syri --prefix $bs -c alns/$bs.paf -r ref.fna.gz -q $seq --lf $bs.syri.log"
	syri --nc 5 -F P --cigar --dir syri --prefix $bs -c alns/$bs.paf -r ref.fna.gz -q $seq --lf $bs.syri.log
done
