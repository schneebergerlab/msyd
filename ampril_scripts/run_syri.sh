#!/bin/sh

for ref in an24 col cvi eri kyo ler sha; do
	for sample in an24 col cvi eri kyo ler sha; do
		echo doing $ref\_$sample
		syri --nc 5 -F B --cigar --dir ~/pansyri/data/ampril --prefix $ref\_$sample -c $ref\_$sample.bam -r $ref.filtered.fa.gz -q $sample.filtered.fa.gz --lf $ref\_$samplesyri.log --samplename $sample
	done
done

