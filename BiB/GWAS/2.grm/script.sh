#!/bin/bash
	
	#PBS -N grm_job
	#PBS -j oe
	#PBS -o bib1-22.out
	#PBS -q himem
	#PBS -l nodes=1:ppn=16
	#PBS -l walltime=48:00:00
	#PBS -t 1-22
	
module add apps/gcta/1.92

gcta --bgen path/to/data_chr0${PBS_ARRAYID}.bgen \
--make-grm-part 5 $${PBS_ARRAYID} --thread-num 5 --maf 0.01 --out ~/grm/bib_cal _cleaned_maf_0_01
