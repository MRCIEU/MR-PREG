#!/bin/bash
#SBATCH --job-name=moba-XXXX-step1-bin-YYYY

#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --mem-per-cpu=10GB


##Set up job environment:
source /cluster/bin/jobsetup
module purge
set -o errexit

module add Regenie/3.1.2-GCC-11.2.0

cd WWWW

regenie --step 1 \
--bed ./data/data-filtered-maf-nomhc-for-regenie \
--covarFile ./data/covs.txt \
--phenoFile ./data/pheno/XXXX_YYYY_v2.txt \
--covarColList genotyping_batch,PC{1:20} \
--catCovarList genotyping_batch \
--maxCatLevels 27 \
--bsize 1000 \
--bt \
--lowmem \
--lowmem-prefix tmp_rg_XXXX_YYYY \
--out ./results/step1_bin_XXXX_YYYY_v2
