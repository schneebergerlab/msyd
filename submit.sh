#!/bin/bash

ref=$1
refpath=$ref/$ref.fna
outpath=${ref}_alns
mkdir -p $outpath
tmpdir="/tmp/lrauschning/"
mkdir -p $tmpdir
for org in albomicans busckii melanogaster sechellia teissieri miranda simulans triauraria subobscura willistoni athabasca lowei pseudoobscura subpulchrella yakuba bifasciata mauritiana santomea suzukii
do
	if [[ "$org" == "$ref" ]]; then
		continue
	fi
	orgpath=${org}/${org}.fna

	bsub -q normal -n4 -R"span[hosts=1] rusage[mem=4000]" -M8000 \
	-oo ${tmpdir}minimap_$ref$org.log -eo ${tmpdir}minimap_$ref$org.err \
	"minimap2 -cx asm20 $refpath $orgpath > $outpath/${ref}_$org.paf"

	# or use asm10?
done
