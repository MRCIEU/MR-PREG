#!/bin/bash
#SBATCH --partition=mrcieu2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=02:00:00
#SBATCH --mem=500M

#SBATCH --array=0-8

source config.sh
cd ${runDir}

module add apps/regenie/3.6

mapfile -t traitnames < outlist.cont.txt
traitname=${traitnames[${SLURM_ARRAY_TASK_ID}]}
prefix=BIB.${traitname}.${personM}.${ethnicitySA}
echo ${traitname}

regenie --step 1 \
--bed ${dat}/step1-saige-both-pcs-sa \
--covarFile ${dat}/pheno.sa.master.txt \
--phenoFile ${dat}/pheno.sa.master.txt  \
--phenoCol ${traitname} \
--covarColList gt_chip_4cat,${covs} \
--catCovarList gt_chip_4cat \
--maxCatLevels 4 \
--bsize 1000 \
--lowmem \
--lowmem-prefix ${lowmem}/${prefix} \
--threads 10 \
--out ${modelsMotherSA}/${prefix}_step1


