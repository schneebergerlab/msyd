#!/bin/bash

bsub -q normal -n4 -R"span[hosts=1] rusage[mem=20000]" -M40000 \
-oo job_all.log -eo job_all.err \
"python ../../pansr/main.py all.tsv 4"
