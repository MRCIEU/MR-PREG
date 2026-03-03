#!/bin/bash
#SBATCH --job-name=moba-dad-step2-cont-rollout

#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10M
#SBATCH --array=1-2

##Set up job environment:
source /cluster/bin/jobsetup
module purge
set -o errexit

## Set working directory
wdir="./MoBaPsychGen_v1/pipeline"
cd ${wdir}/scripts
d=".\/MoBaPsychGen_v1\/pipeline"
echo ${d}

## Define subcohort
i="dad"
echo ${i}

## Index of phenotypes
pIdx=$SLURM_ARRAY_TASK_ID
echo ${pIdx}

## List of phenotypes names
pheno_list=${wdir}/data/pheno/cont_pheno_v2.txt

## Index phenotype 
j=`awk 'NR=='${pIdx}'' ${pheno_list}`
echo ${j}

## Submit GWAS (one job per subcohort & phenotype)
cp step2-cont-pheno_v2.sh step2-${i}-${j}_v2.sh
sed -i 's/WWWW/'${d}'/' step2-${i}-${j}_v2.sh
sed -i 's/XXXX/'${i}'/g' step2-${i}-${j}_v2.sh
sed -i 's/YYYY/'${j}'/g' step2-${i}-${j}_v2.sh
sbatch step2-${i}-${j}_v2.sh
