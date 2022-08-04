#!/bin/bash

bsub -q normal -n4 -R"span[hosts=1] rusage[mem=64000]" -M1280000 \
-oo job_$1.log -eo job_$1.err \
"python ../../pansr/main.py $1 4"
