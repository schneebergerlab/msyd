#!/bin/bash

for n in $(seq 1 20)
do
	bsub -q short -n1 -R"span[hosts=1] rusage[mem=5000]" -M80000 \
	-oo bench_short_$n.log -eo bench_short_$n.err \
	"time pansyri -i ./all.tsv -c 1 --discard"
done
