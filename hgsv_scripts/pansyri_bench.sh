#!/bin/bash

for n in $(seq 1 20)
do
	bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=64000]" -M1280000 \
	-oo bench_$n.log -eo bench_$n.err \
	"pansyri -i ./all.tsv -c 20 --discard"
done
