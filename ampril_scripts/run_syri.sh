#!/bin/bash

for ref in an1 c24 col cvi eri kyo ler sha; do
	for sample in an1 c24 col cvi eri kyo ler sha; do
		if [[ $ref != $sample ]]; then
			echo doing $ref\_$sample
			syri --nc 5 -F B --cigar --dir ~/pansyn/data/ampril --prefix $ref\_$sample -c $ref\_$sample.bam -r $ref.filtered.fa.gz -q $sample.filtered.fa.gz --lf $ref\_$samplesyri.log --samplename $sample
		fi
	done
done

