#!/usr/bin/env bash

#SBATCH --job-name=mrpreg-gwama-child
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=8:00:00
#SBATCH --mem=50G












####################################################################################
##         This is the master script to run the MR-PREG GWAMA for offspring         ##
####################################################################################

# Record some potentially useful details about the job:
echo "Running on host $(hostname)"
echo "Started on $(date)"
echo "Directory is $(pwd)"
echo "Slurm job ID is ${SLURM_JOBID}"
echo "This jobs runs on the following machines:"
echo "${SLURM_JOB_NODELIST}"
printf "\n\n"

echo "Use metal to metanalyse SNP-outcome associations from ALSPAC, BiB, MoBa, FinnGen, UKB and public GWAS"

# Directories
scripts_dir=${MRPREG_sumdat}scripts/GWAS/Metanalysis/offspring/
cd ${scripts_dir}

# List of outcomes for GWAMA
vars_csv=vars_metanalysis_v5_child.csv
awk 'NR>=2' ${vars_csv} | awk -F "\"*,\"*" '{print $2}' > vars_child.txt
vars=$(cat vars_child.txt)
echo ${vars}

# Write METAL scripts for each outcome
module add languages/R/4.3.3 
R CMD BATCH 1_write_metal_scripts_child.R 1_write_metal_scripts_child.Rout

# Run METAL
module add metal/2020-05-05 
declare -a phenos=${vars}
for pheno in ${phenos[@]}; do 
  ## Print phenotype
   echo ${pheno} 
  ## Run METAL scripts
   metal ./Linux/metal_${pheno}.txt
done

# Format METAL output
R CMD BATCH 2_format_metal_output_child.R 2_format_metal_output_child.Rout

# Print ending time:
printf "\n\n"
echo "Ended on: $(date)"