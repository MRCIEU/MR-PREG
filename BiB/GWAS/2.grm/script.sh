#!/bin/bash

#SBATCH --job-name=make-grm
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=24:0:0
#SBATCH --mem=50gb
#SBATCH --array=1-22

module add apps/gcta/1.92

gcta --bgen path/to/data_chr${SLURM_ARRAY_TASK_ID}.bgen \
--make-grm-part 50 ${SLURM_ARRAY_TASK_ID} \
--thread-num 5 \
--maf 0.01 \
--out path/to/bib_cal_cleaned_maf_0_01


