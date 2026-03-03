#!/bin/bash
#SBATCH --job-name=step2_cont_[1-5]
#SBATCH --partition=mrcieu2

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=4GB
#SBATCH --time=24:00:00
#SBATCH --array=1-5

## set up job
echo "Running on ${HOSTNAME}"
module load apps/regenie/3.6
ukbdir=$UKB_DIR/released/2018-09-18/data/dosage_bgen
scratchdir=$MRPREG_sumdat

## list of phenotype names
pheno_list=${scratchdir}/data/GWAS/UKB/mothers/cont_pheno.txt

## array index
pIdx=$SLURM_ARRAY_TASK_ID
echo ${pIdx}

## index phenotype
pheno=`awk 'NR=='${pIdx}'' ${pheno_list} `
echo ${pheno}

## Run regenie step 2 (per pheno, per chr)
for chr in {1..22}; do
### zero pad chromosomes <10 to match bgen file names
if [ ${chr} -lt 10 ]; then
chr="0"${chr}
fi
### Step 2
regenie \
  --step 2 \
  --bgen ${ukbdir}/data.chr${chr}.bgen \
  --covarFile ${scratchdir}/data/GWAS/UKB/mothers/covs.txt \
  --phenoFile ${scratchdir}/data/GWAS/UKB/mothers/pheno_files/cont_${pheno}.txt \
  --sample ${ukbdir}/data.chr1-22.sample \
  --bsize 200 \
  --qt \
  --maxCatLevels 27 \
  --threads 10 \
  --pred ${scratchdir}/results/GWAS/UKB/mothers/step1_cont_${pheno}_pred.list \
  --out ${scratchdir}/results/GWAS/UKB/mothers/step2_cont_${pheno}_chr${chr}
done