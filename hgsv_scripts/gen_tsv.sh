#!/bin/sh

echo "#name aln 	syri	vcf" > all.tsv
for f in syri/*syri.out
do
	bs=$(basename $f syri.out)
	# shorten the name even further to make it readable
	name="$(echo $bs | grep -o -P '(NA|HG)\d+').$(echo $bs | grep -o -P 'h[1|0]')"
	echo "$name alns_filtered/$bs.paf	syri/${bs}syri.out syr/${bs}syri.vcf" >> all.tsv
done
