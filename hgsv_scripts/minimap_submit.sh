#!/bin/sh

ref=$1
outpath="alns_filtered"

for seq in scaffolds_filtered/*.fna.bgz
do
	bs=$(basename $seq)
	bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=20000]" -M30000 \
		-oo aln_$bs.log -eo aln_$bs.err \
		"minimap2 -cx asm5 --eqx ref-nomtsex.fna.bgz $seq > $outpath/$bs.paf"
done
