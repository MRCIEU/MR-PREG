## GRAMMAR-Gamma is an approximate method that gives O(MN) time of SNP testing instead of O(MN2) 
##(M is the number of SNPs, N - that of individuals)

#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

module load apps/gcta/1.93.2-beta

gcta64 --bgen /bp1/mrcieu1/data/bib/genetic/variants/arrays/imputed/uk10k_1000genomes/released/2017-01-27/data/bgen/tom/data_chr${SLURM_ARRAY_TASK_ID}.bgen --grm-sparse sp_grm --fastGWA-lr --pheno /user/home/rt20627/pheno.txt --out testrun
## we got an error of not finding the pheno file and we suspect that we need to specigy the grms (--grm-sparse "sp_grm" to be changed). we tried by reordering the columns of the pheno file and that did not work

#--fastGWA-mlm --data --est-vg --h2-limit 1.6 --save-fastGWA-mlm-residual --seed 123


#run in the same sample
--fastGWA-lr linear regression
--fastGWA-mlm linear regression
#Generate histogram of pvals to see if the distribution is like Tom’s (normal) or uniform (as expected) – check with GWAS inspector package 
