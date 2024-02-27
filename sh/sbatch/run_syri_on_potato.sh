#!/bin/bash
#SBATCH --array=0-11
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5000mb
#SBATCH --time=12:00:00
#SBATCH -J QUickScript

#module load samtools

cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/msyd/results/potato_hap/
indir=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/msyd/data/potato_haps/
cd $cwd

chars=({A..J} O)
for r in 1 2 3 4; do
  for q in 1 2 3 4; do
    echo ${chars[${SLURM_ARRAY_TASK_ID}]}
    srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=5000 \
      /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam \
      -n 10 -p O${r}_vs_${chars[${SLURM_ARRAY_TASK_ID}]}${q} \
      ${indir}/O_hap${r}_genome.fasta ${indir}/${chars[${SLURM_ARRAY_TASK_ID}]}_hap${q}_genome.fasta &
  done
done
wait

#for i in 5 6 7 8; do
#
#	echo ${SLURM_CPUS_PER_TASK}
#	echo ${chars[${SLURM_ARRAY_TASK_ID}]}
#	srun --exclusive --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem-per-cpu=5000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam \
#	-n 10 -p dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i} \
#	DM_chr02.fa ${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i}.fa &
##srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=1000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam -n 10 -p chr02_dm_A_hap${i} DM_chr02.fa chr02_hap${i}.fa
#done
#wait

#<<interactive
#for c in {A..J}; do
#for c in {H..H}; do
#  for i in 5 6 7 8; do
#  #  srun minimap2 -ax asm5 -t 10 --eqx DM_chr02.fa ${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i}.fa \
#  #  | samtools sort -@ 10 -O BAM -o dm_${chars[${SLURM_ARRAY_TASK_ID}]}_chr02_hap${i}.bam -
#  #  srun samtools
#    srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=5000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam \
#    -n 10 -p dm_${c}_chr02_hap${i} \
#    DM_chr02.fa ${c}_chr02_hap${i}.fa &
#  #srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=1000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam -n 10 -p chr02_dm_A_hap${i} DM_chr02.fa chr02_hap${i}.fa
#  done
#  wait
#done
##wait
#interactive




#module load samtools
##chars=({O..O})
##for i in 5 6 7 8; do
#srun /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools \
#  runsyri -alignment bam \
#	-n 10 -p dm_O_chr02_hap${SLURM_ARRAY_TASK_ID} \
#	DM_chr02.fa O_chr02_hap${SLURM_ARRAY_TASK_ID}.fa
##srun --ntasks 1 --cpus-per-task ${SLURM_CPUS_PER_TASK} --mem-per-cpu=1000 /dss/dsslegfs01/pn29fi/pn29fi-dss-0003/software/bin_manish/anaconda3/envs/mgpy3.8/bin/hometools runsyri -alignment bam -n 10 -p chr02_dm_A_hap${i} DM_chr02.fa chr02_hap${i}.fa
##done
##wait
#

