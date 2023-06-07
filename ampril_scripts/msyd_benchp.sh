#!/bin/bash
# benchmark the parallelizability of the algorithm

for c in $(seq 1 10)
do
	for n in $(seq 1 5)
	do
		bsub -q multicore20 -n$c -R"span[hosts=1] rusage[mem=8000]" -M320000 \
		-oo benchp_$c_$n.log -eo benchp_$c_$n.err \
		"msyd call -i ./full.tsv -c $c -o /dev/null"
	done
done
