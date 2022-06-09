#!/bin/sh
tmp="/tmp/"
out="scaffolds_filtered"
file=$(basename $1 .bgz)

# determine chr list from reference
# mt-removed ref stored in basedir as ref-nomt.fna.bgz
chrs=$(zgrep -e ">NC" grch38.v14_genomic.fna.gz | sed "s/>//" | cut -f 1 -d " " | head -n 24 | sed -z "s/\n/\ /g")

zcat scaffolds/$file.bgz | sed "s/_RagTag$//" > $tmp/$file.tmp

hometools getchr --chrs $chrs -o scaffolds_filtered/$file $tmp/$file.tmp
rm $tmp/$file.tmp
bgzip --threads 4 scaffolds_filtered/$file
mv scaffolds_filtered/$file.gz scaffolds_filtered/$file.bgz
