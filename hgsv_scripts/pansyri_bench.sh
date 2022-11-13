#!/bin/bash

for n in $(seq 1 20)
do
	bsub -q multicore20 -n10 -R"span[hosts=1] rusage[mem=128000]" -M1280000 \
	-oo bench_new_$n.log -eo bench_new_$n.err \
	"pansyri -i ./haplos.tsv -c 10 --discard"
done
