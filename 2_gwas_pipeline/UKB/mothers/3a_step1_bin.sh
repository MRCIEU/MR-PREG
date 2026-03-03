#!/bin/bash
#SBATCH --job-name=step1_bin_[1-42]
#SBATCH --partition=mrcieu

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=10GB
#SBATCH --time=36:00:00
#SBATCH --array=1-42

## set up job
echo "Running on ${HOSTNAME}"

module load apps/regenie/3.6

ukbdir=$UKB_DIR/released/2017-07-04/data/derived/merged_chr1-22
scratchdir=$MRPREG_sumdat

## index of phenotypes
pIdx=$SLURM_ARRAY_TASK_ID
echo ${pIdx}

## list of phenotype names
pheno_list=${scratchdir}/data/GWAS/UKB/mothers/bin_pheno.txt

## index phenotype
j=`awk 'NR=='${pIdx}'' ${pheno_list} `
echo ${j}

regenie \
  --step 1 \
  --bed ${ukbdir}/chr1-22_merged \
  --extract ${scratchdir}/results/GWAS/UKB/mothers/qc-pass.snplist \
  --phenoFile ${scratchdir}/data/GWAS/UKB/mothers/pheno_files/bin_${j}.txt \
  --covarFile ${scratchdir}/data/GWAS/UKB/mothers/covs.txt \
  --keep ${scratchdir}/results/GWAS/UKB/mothers/qc-pass.id \
  --bt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_rg_${j} \
  --threads 10 \
  --out ${scratchdir}/results/GWAS/UKB/mothers/step1_bin_${j}