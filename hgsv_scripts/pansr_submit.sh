#!/bin/bash

bsub -q normal -n4 -R"span[hosts=1] rusage[mem=64000]" -M1280000 \
-oo job_all.log -eo job_all.err \
"python ../../pansr/main.py $1 4"
