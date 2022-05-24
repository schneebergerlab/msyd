#!/bin/sh

# Ensure PATH contains /netscratch/dep_mercier/grp_schneeberger/software/RagTag-2.1.0/
tmpdir="/tmp/lrauschning/"
for i in $1/*.fasta.bgz
do
	bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=32000]" -M32000 \
		-oo ragtag_$(basename $i).log -eo ragtag_$(basename $i).err \
		"ragtag.py scaffold grch38.v14_genomic.fna.gz $i -t 20 -o scaffolds/$(basename $i)"# -w"
done
