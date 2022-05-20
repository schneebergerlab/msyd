#!/bin/sh

# set up to work with the files in data/drosophila

cd $1
for f in *.fna; do
    base=$(basename "$f" '.fna') # gives '25' from '25.conf'
    sed -i "1 s/.*/>$1_$base/" "$f"
done

mv unplaced.scaf.fna unplaced.scaf
cat *.fna > $1.fna
chmod -x $1.fna
