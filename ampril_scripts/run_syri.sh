#!/bin/sh

for ref in an24 col cvi eri kyo ler sha; do
	echo doing $ref
	for sample in an24 col cvi eri kyo ler sha; do
		syri --nc 5 -F B --cigar --dir ~/pansyri/data/ampril --prefix $ref_$sample -c $ref_$sample.bam -r $ref.filtered.fa.gz -q $sample.filtered.fa.gz --lf $ref_$samplesyri.log --samplename $sample
	done
done

