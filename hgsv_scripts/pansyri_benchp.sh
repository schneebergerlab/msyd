#!/bin/bash
# benchmark the parallelizability of the algorithm

for c in $(seq 20 4 40)
do
	for n in $(seq 1 5)
	do
		bsub -q multicore40 -n$c -R"span[hosts=1] rusage[mem=128000]" -M1280000 \
		-oo benchp_$c_$n.log -eo benchp_$c_$n.err \
		"pansyri -i ./haplos.tsv -c $c --discard"
	done
done
