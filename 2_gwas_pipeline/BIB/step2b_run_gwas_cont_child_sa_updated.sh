#!/bin/bash
#SBATCH --partition=mrcieu2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=04:00:00
#SBATCH --mem=1G

#SBATCH --array=0-8

source config.sh
cd ${runDir}

module add apps/regenie/3.6

mapfile -t traitnames < outlist.cont.txt
traitname=${traitnames[${SLURM_ARRAY_TASK_ID}]}
prefix=BIB.${traitname}.${personC}.${ethnicitySA}
echo ${traitname}

for i in {1..22}; do
  regenie --step 2 \
  --bgen ${datSA}filtered_R203_chr${i}_zlib.bgen \
  --sample ${datSA}filtered_R203.sample \
  --bgi ${datSA}filtered_R203_chr${i}_zlib.bgen.bgi \
  --covarFile ${dat}/pheno.sa.master.children.txt \
  --phenoFile ${dat}/pheno.sa.master.children.txt \
  --phenoCol ${traitname} \
  --covarColList gt_chip_4cat,${covs} \
  --catCovarList gt_chip_4cat \
  --maxCatLevels 4 \
  --bsize 200 \
  --minMAC 20 \
  --threads 10 \
  --pred  ${modelsOffSA}/${prefix}_step1_pred.list \
  --out ${resOffSA}/${prefix}_chr${i}
done
