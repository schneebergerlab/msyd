#!/bin/bash

for n in $(seq 1 20)
do
	bsub -q short -n4 -R"span[hosts=1] rusage[mem=20000]" -M320000 \
	-oo bench_new_$n.log -eo bench_new_$n.err \
	"time pansyri -i ./full.tsv -c 4 --discard"
done
