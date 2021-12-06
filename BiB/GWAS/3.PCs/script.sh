
#!/bin/bash

#SBATCH --job-name=pcs_job
#SBATCH --partition=hmem
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22


#calculate pcs for unrelated individuals, then project all individuals including related individuals onto these pcs
	

#define variables 
gen_files=/path/to/vcf/
hardcalled_files=/path/to/hardcalled/chr${SLURM_ARRAY_TASK_ID}

out_files=path/to/out/
set_file=price_2008_long_range_ld_snps.txt 
#check if SLURM_ARRAY_TASK_ID<10
if [ ${SLURM_ARRAY_TASK_ID} -gt 9 ]
then 
	gen_file=${gen_files}data_chr${SLURM_ARRAY_TASK_ID}.vcf.gz
else 
	gen_file=${gen_files}data_chr0${SLURM_ARRAY_TASK_ID}.vcf.gz
fi

echo $gen_file

hardcalled_file=${hardcalled_files}/chr${SLURM_ARRAY_TASK_ID}
out_file=${out_files}chr${SLURM_ARRAY_TASK_ID}
plink=path/to/plink
plink2=path/to/plink2
#flashpca=


#to write a set of plink bed/bim/fam files, which contain “hard called” genotype at imputed snps (i.e. integer values representing the most likely imputed genotype).
time ${plink2} --vcf $gen_file dosage=GP --mach-r2-filter 0.82 --hard-call-threshold 0.1 --import-dosage-certainty 0.9 --make-bed --out $hardcalled_file
 

#make .set file with a list of variants in long range ld (from Price, Alkes L., et al. "Long-range LD can confound genome scans in admixed populations." The American Journal of Human Genetics 83.1 (2008): 132-135.)
time ${plink} --bfile $hardcalled_file \
 --make-set $set_file \
 --write-set \
 --threads 20 \
 --out ${out_files}long_range_ld${SLURM_ARRAY_TASK_ID}
 
 # remove related individuals
${plink2} \
--bfile ${hardcalled_file} \
--autosome \
--threads 20 \
--king-cutoff 0.04419417 \
--make-bed \
--out ${out_file}_unrelated

# prune snps (flashPCA defaults (https://github.com/gabraham/flashpca), having removed relateds, and remove long range ld
${plink2} \
--bfile ${out_file}_unrelated \
--autosome \
--maf 0.01 \
--indep-pairwise 1000 50 0.05 \
--threads 20 \
--exclude ${out_files}long_range_ld${SLURM_ARRAY_TASK_ID}.set \
--out ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned

# write plink binary files to calculate pcs from
${plink2} \
--bfile ${out_file}_unrelated \
--make-bed \
--extract ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned.prune.in \
--out ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned

# intermediate files can be deleted
rm ${out_file}_unrelated.bed ${out_file}_unrelated.bim ${out_file}_unrelated.fam

# calculate pcs using flashPCA (need to check what command to use to get 100 pcs; see https://github.com/gabraham/flashpca)
time ./flashpca --bfile ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned -d 100 \
--outload ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_loadings.txt \
--outmeansd ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_meansd.txt 


#####TBC#####
