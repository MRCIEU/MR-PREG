#!/bin/bash
#SBATCH --partition=mrcieu2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00


cd ${MRPREG_sumdat}/

# Raw pheno file
#rsync ${RDSF_IEU2}/data/UKB/outcomes/Phenotype_MRPREG.txt ./data/GWAS/UKB/mothers/
#rsync ${RDSF_IEU2}/data/UKB/outcomes/phenotypic_data_maternity.txt ./data/GWAS/UKB/mothers/

# WD
cd ./scripts/GWAS/UKB/mothers/

# Load R
module add languages/R/4.3.3

# Prepare pheno file
R CMD BATCH 1b_separate_phenotypes.R 1b_separate_phenotypes.Rout

# Prepare covariables file
R CMD BATCH 1c_combine_covariates.R 1c_combine_covariates.Rout