#!/bin/sh

ref=$1
outpath="alns"

for seq in scaffolds/*.fna.gz
do
	bs=$(basename $seq)
	bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=20000]" -M30000 \
		-oo aln_$bs.log -eo aln_$bs.err \
		"minimap2 -cx asm5 --eqx grch38.v14_genomic.fna.gz $seq > $outpath/$bs.paf"
done
