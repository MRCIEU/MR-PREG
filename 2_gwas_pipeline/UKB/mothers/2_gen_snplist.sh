#!/bin/bash
#SBATCH --partition=mrcieu

#SBATCH --job-name=plink_step0
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

echo "Running on ${HOSTNAME}"

cd $UKB_DIR

#module load apps/plink/2.3
module add apps/plink2/2.00a6LM

plink2 \
  --bfile ./released/2017-07-04/data/derived/merged_chr1-22/chr1-22_merged \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --keep-females \
  --write-snplist --write-samples --no-id-header \
  --out $MRPREG_sumdat/results/GWAS/UKB/mothers/qc-pass
