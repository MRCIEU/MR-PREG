#!/bin/bash

#SBATCH --job-name=pcapred_job
#SBATCH --partition=mrcieu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=300000M

eval $R
time Rscript pcapred.R
