#!/bin/sh

echo "# aln 	syri" > all.tsv
for f in syri/*syri.out
do
	bs=$(basename $f syri.out)
	echo "alns/$bs.paf	syri/${bs}syri.out" >> all.tsv
done
