#!/bin/bash

out=$(basename $1 .pff)

bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=64000]" -M1280000 \
-oo job_$1.log -eo job_$1.err \
"pansyn call -i $1 -c 20 -o $out.pff -m $out.vcf -r ref-nomtsex.fna.bgz"
