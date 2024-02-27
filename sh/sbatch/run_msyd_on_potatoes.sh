#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000mb
#SBATCH --time=12:00:00
#SBATCH -J QUickScript

module load samtools


#for i in 5 6 7 8; do
  srun  /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/msyd call --realign -i genomes.csv -o chr02_all_hap.pff -r DM_chr02.fa
#srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem=10000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam -n 10 -p chr02_dm_A_hap${i} DM_chr02.fa chr02_hap${i}.fa
#done
#wait
