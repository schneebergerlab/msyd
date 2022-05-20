#!/bin/bash

ref=$1
refpath=sequences/$ref/$ref.fna
outpath=${ref}_alns
mkdir -p $outpath
tmpdir="/tmp/lrauschning/"
mkdir -p $tmpdir

# all flies:
all="albomicans busckii melanogaster sechellia teissieri miranda simulans triauraria subobscura willistoni athabasca lowei pseudoobscura subpulchrella yakuba bifasciata mauritiana santomea suzukii"
# flies that have the same chromosome structure as melanogaster/simulans:
sim="busckii melanogaster sechellia teissieri simulans triauraria subpulchrella yakuba mauritiana santomea"
for org in $sim
do
	if [[ "$org" == "$ref" ]]; then
		continue
	fi
	orgpath=sequences/${org}/${org}.fna

	bsub -q normal -n5 -R"span[hosts=1] rusage[mem=4000]" -M8000 \
	-oo ${tmpdir}minimap_$ref$org.log -eo ${tmpdir}minimap_$ref$org.err \
	"minimap2 -cx asm20 --eqx $refpath $orgpath > $outpath/${ref}_$org.paf"

	# or use asm10?
done
