#!/bin/bash
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=02:00:00
#SBATCH --mem=4G

#SBATCH --array=0-40

source config.sh
cd ${runDir}

module add apps/regenie/3.6

mapfile -t traitnames < outlist.bin.txt
traitname=${traitnames[${SLURM_ARRAY_TASK_ID}]}
prefix=BIB.${traitname}.${personC}.${ethnicitySA}
echo ${traitname}

regenie --step 1 \
--bed ${dat}/step1-saige-both-pcs-sa \
--covarFile ${dat}/pheno.sa.master.children.txt \
--phenoFile ${dat}/pheno.sa.master.children.txt  \
--phenoCol ${traitname} \
--covarColList gt_chip_4cat,${covs} \
--catCovarList gt_chip_4cat \
--maxCatLevels 4 \
--bsize 1000 \
--bt \
--lowmem \
--lowmem-prefix ${lowmem}/${prefix} \
--threads 10 \
--out ${modelsOffSA}/${prefix}_step1


