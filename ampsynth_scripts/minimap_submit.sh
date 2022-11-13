#!/bin/sh

ref="./ref.fna.gz"
outpath="alns"

for seq in seqs/*.fna
do
	bs=$(basename -s .fna.gz $seq)
	bsub -q multicore20 -n5 -R"span[hosts=1] rusage[mem=10000]" -M20000 \
		-oo aln_$bs.log -eo aln_$bs.err \
		"minimap2 -cx asm5 --eqx ref.fna.gz $seq > $outpath/$bs.paf"
done
