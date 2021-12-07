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

gcta --bgen /mnt/storage/private/mrcieu/data/bib/genetic/variants/arrays/imputed/uk10k_1000genomes/released/2017-01-27/data/bgen/tom/data_chr${SLURM_ARRAY_TASK_ID}.bgen \
--make-grm-part 5 1 \
--thread-num 5 \
--maf 0.01 \
--out path/to/bib_cal_cleaned_maf_0_01


