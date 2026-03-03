#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=50gb









## Load config.sh file for paths 

source config.sh 

## White European 

plink2 --bfile ${bib_dg} \
--maf 0.01 \
--hwe 0.000001 \
--geno 0.03 \
--keep ${datE}/filtered_R203.sample \
--autosome \
--write-snplist \
--make-bed \
--out ${datE}/step1-saige-both-pcs-white

#awk '{$2 = $1 " " $2}1' step1-saige-both-pcs-white.fam > step1-saige-white-2.fam

## South Asian 

plink2 --bfile ${bib_dg} \
--maf 0.01 \
--hwe 0.000001 \
--geno 0.03 \
--keep ${datSA}/filtered_R203.sample \
--autosome \
--write-snplist \
--make-bed \
--out ${datSA}/step1-saige-both-pcs-sa

#awk '{$2 = $1 " " $2}1' step1-saige-both-pcs-sa.fam > step1-saige-sa-2.fam 
