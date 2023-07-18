#!/bin/sh

echo "#name aln 	syri	vcf	sequence" > all.tsv
for f in syri/*syri.out
do
	bs=$(basename $f syri.out)
	# shorten the name even further to make it readable
	name="$(echo $bs | grep -o -P '(NA|HG)\d+').$(echo $bs | grep -o -P 'h[2|1|0]').$(echo $bs | grep -o -P 'flye')"
	echo "$name	alns_filtered/$bs.paf	syri/${bs}syri.out	syri/${bs}syri.vcf	scaffolds_filtered/$bs" >> all.tsv
done
