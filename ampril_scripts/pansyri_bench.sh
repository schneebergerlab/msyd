#!/bin/bash

for n in $(seq 1 20)
do
	bsub -q multicore20 -n8 -R"span[hosts=1] rusage[mem=8000]" -M320000 \
	-oo bench_$n.log -eo bench_$n.err \
	"pansyri -i ./full.tsv -c 8 --discard"
done
