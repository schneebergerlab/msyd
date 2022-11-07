#!/bin/bash

bsub -q multicore20 -n20 -R"span[hosts=1] rusage[mem=64000]" -M1280000 \
-oo job_$1.log -eo job_$1.err \
"pansyri -i $1 -c 20 --sp"
