#!/bin/bash

#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mem=10gb
#SBATCH --array=1-11


#---------------------------------------------

##
## Record some useful details about the job:
echo Running on host `hostname`
echo Time is `date`
echo SLURMjob ID is  $SLURM_JOB_ID
echo "MR-PREG outcomes GWAS - Mothers"

##
## Scripts directory
cd ${MRPREG_sumdat}scripts/GWAS/ALSPAC/mothers/


##
## Start index of phenotypes
pIdx=${SLURM_ARRAY_TASK_ID}
echo ${pIdx}


##
## List of traits (ie MR-PREG outcomes) 
outlist=${MRPREG_sumdat}data/GWAS/ALSPAC/mothers/cont_pheno.txt


##
## Index trait 
phen=`awk 'NR=='${pIdx}'' ${outlist}`
echo ${phen}


##
## Run GWAS

# Imputed genotype dataset (.bgen) released 04.05.2017
gendir="${ALSPAC_HRC}/bgen"
# Sample file
sampfile="${MRPREG_sumdat}data/GWAS/ALSPAC/mothers/alspac_mums.sample"
# Phenotype dataset
phenofile="${MRPREG_sumdat}data/GWAS/ALSPAC/mothers/alspac_mum.phen"
# Covariables dataset
covarfile="${MRPREG_sumdat}data/GWAS/ALSPAC/mothers/alspac_mums.pcs"
# List of people that have withdrawn consent (last updated 03/10/2023) - Already excluded in the phenotype dataset 
# List of related individuals to exclude
related="${ALSPAC_GENO}/derived/unrelated_ids/mothers/exclusion_list.txt"     
# Output file
outdir="${MRPREG_sumdat}results/GWAS/ALSPAC/mothers/raw"


for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22; do
outfile=mums_plink_${phen}_chr$i
#${MRPREG_sumdat}tmp/plink2 \
${MRPREG_sumdat}/software/plink2_linux_avx2_20241222 \
--bgen ${gendir}/data_$i.bgen ref-unknown \
--sample ${sampfile} \
--pheno ${phenofile} \
--covar ${covarfile} \
--pheno-name ${phen} \
--remove ${related} \
--threads 8 \
--glm hide-covar \
--debug \
--out ${outdir}/${outfile}
done

