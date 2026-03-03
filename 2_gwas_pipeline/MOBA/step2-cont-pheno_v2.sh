#!/bin/bash
#SBATCH --job-name=moba-XXXX-step2-cont-YYYY

#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --mem-per-cpu=10GB

##Set up job environment:
source /cluster/bin/jobsetup
module purge
set -o errexit

module add Regenie/3.1.2-GCC-11.2.0

cd WWWW

regenie --step 2 \
--bed ./MoBaPsychGen_v1-ec-eur-batch-basic-qc \
--covarFile ./pipeline/data/covs.txt \
--phenoFile ./pipeline/data/pheno/XXXX_YYYY_v2.txt \
--covarColList genotyping_batch,PC{1:20} \
--catCovarList genotyping_batch \
--maxCatLevels 27 \
--bsize 200 \
--pred ./pipeline/results/step1_cont_XXXX_YYYY_v2_pred.list \
--out ./pipeline/results/step2_cont_XXXX_YYYY_v2
