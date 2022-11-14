#!/bin/bash

for n in $(seq 1 20)
do
	bsub -q multicore20 -n4 -R"span[hosts=1] rusage[mem=5000]" -M80000 \
	-oo job_$n.log -eo job_$n.err \
	"pansyri -i ./all.tsv -c 4 --discard"
done
