#!/bin/bash

#PBS -N bgen_job
#PBS -j oe
#PBS -o bib1-22.out
#PBS -q himem
#PBS -l nodes=1:ppn=16
#PBS -l walltime=48:00:00
#PBS -t 1-22

module add apps/qctool-2.1

echo "Current working directory is `pwd`" or $PBS_O_WORKDIR
echo "Running on `hostname`" 
echo job ID is $PBS_JOBID 
echo "But output is bib1-22.out"
echo This jobs runs on nodes `cat $PBS_NODEFILE | uniq`
echo BC3 has 16 high-memory nodes, each of which has 256GB of RAM. We request 16 processors in a single compute node. 
echo count the number of processors available: numprocs=`wc $PBS_NODEFILE | awk '{print $1}'`
echo "job started at `date`"

qctool -g path/to/data_chr0${PBS_ARRAYID}.vcf.gz -og path/to/data_chr0${PBS_ARRAYID}_bib.bgen -os path/to/sampleid_bib0${PBS_ARRAYID}.sample -sex-column sex
bgenix -index -g path/to/data_chr0${PBS_ARRAYID}_bib.bgen

echo "job finished at `date`"
