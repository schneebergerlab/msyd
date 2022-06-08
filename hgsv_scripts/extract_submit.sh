#!/bin/bash

for seq in scaffolds/*.fna.bgz
do
	bs=$(basename $seq)
	bsub -q normal -n4 -R"span[hosts=1] rusage[mem=8000]" -M10000 \
	-oo job_extract_${bs}.log -eo job_extract_${bs}.err \
	"~/pansr/pansr/hgsv_scripts/extract_chrs.sh $seq"
done
