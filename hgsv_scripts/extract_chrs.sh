#!/bin/sh
tmp="/tmp/"
out="scaffolds_filtered"
file=$(basename $1)

# determine chr list from reference
chrs=$(zgrep -e ">NC" grch38.v14_genomic.fna.gz | sed "s/>//" | cut -f 1 -d " " | sed -z "s/\n/\ /g")

zcat scaffolds/$file | sed "s/_RagTag$//" > $tmp/$file.tmp

hometools getchr --chrs $chrs -o scaffolds_filtered/$file $tmp/$file.tmp
rm $tmp/$file.tmp
bgzip --threads 4 scaffolds_filtered/$file

