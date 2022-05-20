#!/bin/bash

ref=$1
refpath=sequences/$ref/$ref.fna
outpath=${ref}_alns
alnpath=$outpath
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
	aln=$alnpath/${ref}_$org.paf

	bsub -q multicore20 -n6 -R"span[hosts=1] rusage[mem=6000]" -M10000 \
	-oo ${tmpdir}syri_$ref$org.log -eo ${tmpdir}syri_$ref$org.err \
	"syri --nc 6 -F P --cigar --dir $outpath --prefix ${ref}_$org -c $aln -r $refpath -q $orgpath"

	# or use asm10?
done
