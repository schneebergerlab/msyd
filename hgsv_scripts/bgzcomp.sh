#!/bin/sh

tmpdir="/tmp/lrauschning/"
cd scaffolds
for dir in $(ls -d */)
do
	echo copying from $dir to $(basename $dir).fna.bgz
	bsub -q ioheavy -n3 -R"span[hosts=1] rusage[mem=1000]" -M2000 \
		-oo comp_$(basename $dir).log -eo comp_$(basename $dir).err \
		"bgzip -c $dir/ragtag.scaffolds.fasta > $(basename $dir).fna.bgz"
done
