#!/bin/bash

#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=8:00:00
#SBATCH --mem=10gb


##
## Record some useful details about the job:
echo Running on host `hostname`
echo Time is `date`
echo SLURMjob ID is  $SLURM_JOB_ID
echo This jobs runs on the following machines:
echo `cat $SLURM_JOB_NODELIST | uniq`
echo "MR-PREG outcomes GWAS - Mothers - extract EAF"


##
## Load Plink2
#module add plink2/2.00a4.3-openblas


##
## Run GWAS: outcome YYY - all chromosomes

# Imputed genotype dataset (.bgen) released 04.05.2017
gendir="${ALSPAC_HRC}/bgen"
# Sample file
sampfile="${MRPREG_sumdat}data/GWAS/ALSPAC/mothers/alspac_mums.sample"
# List of related individuals to exclude
unrelated="${ALSPAC_GENO}/derived/unrelated_ids/mothers/inclusion_list.txt"     
# Output file
outdir="${MRPREG_sumdat}results/GWAS/ALSPAC/mothers/raw"


for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22; do
${MRPREG_sumdat}tmp/plink2 --bgen ${gendir}/data_$i.bgen ref-unknown --sample ${sampfile} --freq --keep ${unrelated} --out ${outdir}/mums_freq_chr$i
done
