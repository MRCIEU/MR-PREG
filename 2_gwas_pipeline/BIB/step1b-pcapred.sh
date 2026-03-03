#!/bin/bash
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=50gb


source config.sh

Rscript step1b-pcapred.R