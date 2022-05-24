#!/bin/sh

tmpdir="/tmp/lrauschning/"
for i in $1/*.fasta.gz
do
	bsub -q ioheavy -n3 -R"span[hosts=1] rusage[mem=1000]" -M2000 \
		-oo conv_$(basename $i .gz).log -eo conv_$(basename $i).err \
		"pigz -kdc $i | bgzip -c > $(basename $i .gz)"
done
