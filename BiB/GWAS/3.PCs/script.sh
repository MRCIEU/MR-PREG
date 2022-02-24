#Modified by Marwa Al Arab, based on script from Tom Bond.
#!/bin/bash

#SBATCH --job-name=pcs_job
#SBATCH --partition=mrcieu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4


#array=1-22


#calculate pcs for unrelated individuals, then project all individuals including related individuals onto these pcs
#run from /mnt/storage/private/mrcieu/data/bib/genetic/variants/arrays/imputed/uk10k_1000genomes/released/2017-01-27/scripts/pcs	

#define variables 

plink1='module load apps/plink/1.90'
plink2='module load apps/plink/2.00'


#gen_files=/projects/MRC-IEU/research/data/bib/genetic/variants/arrays/directly_genotyped/released/2020-03-18/data/data
gen_files=/user/home/cu20932/bib/directly_genotyped/2020-03-18/data/data

#hardcalled_files = ../../data/pcs/hardcalled/chr${SLURM_ARRAY_TASK_ID}
#projects/MRC-IEU/research/data/bib/genetic/variants/arrays/directly_genotyped/released/2020-03-18/data/data

out_files=../../data/pcs/output/
#out_path=
set_file=price_2008_long_range_ld_snps.txt 
white_mums_ids=included/white.mums.txt
pakis_mums_ids=included/pakis.mums.txt
white_child_ids=included/white.child.txt

#check if SLURM_ARRAY_TASK_ID<10
# if [ ${SLURM_ARRAY_TASK_ID} -gt 9 ]
# then 
	# gen_file=${gen_files}data_chr${SLURM_ARRAY_TASK_ID}.vcf.gz
# else 
	# gen_file=${gen_files}data_chr0${SLURM_ARRAY_TASK_ID}.vcf.gz
# fi

# echo $gen_file

#hardcalled_file=${hardcalled_files}/chr${SLURM_ARRAY_TASK_ID}
#out_file=${out_files}white_mums
#out_file=${out_files}pakis_mums
out_file=${out_files}white_child





#to write a set of plink bed/bim/fam files, which contain “hard called” genotype at imputed snps (i.e. integer values representing the most likely imputed genotype).
#module load apps/plink/2.00
#time plink --vcf $gen_file dosage=GP --mach-r2-filter 0.82 --hard-call-threshold 0.1 --import-dosage-certainty 0.9 --make-bed --out $hardcalled_file
 
eval $plink1
#module load apps/plink/1.90
#make .set file with a list of variants in long range ld (from Price, Alkes L., et al. "Long-range LD can confound genome scans in admixed populations." The American Journal of Human Genetics 83.1 (2008): 132-135.)
 time plink --bfile ${gen_files} \
 --keep ${white_child_ids} \
 --make-set ${set_file} \
 --write-set \
 --threads 20 \
 --out ${out_file}_long_range_ld
 
#remove related individuals
module load apps/plink/2.00
eval $plink2
time plink \
--bfile ${gen_files} \
--autosome \
--threads 20 \
--keep ${white_child_ids} \
--king-cutoff 0.04419417 \
--make-bed \
--out ${out_file}_unrelated


#prune snps (flashPCA defaults (https://github.com/gabraham/flashpca), having removed relateds, and remove long range ld
eval $plink2
plink \
--bfile ${out_file}_unrelated \
--autosome \
--maf 0.01 \
--indep-pairwise 1000 50 0.05 \
--threads 20 \
--exclude ${out_file}_long_range_ld.set \
--out ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned

#write plink binary files to calculate pcs from
eval $plink1
plink \
--bfile ${out_file}_unrelated \
--make-bed \
--extract ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned.prune.in \
--out ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned

#intermediate files can be deleted
 rm ${out_file}_unrelated.bed ${out_file}_unrelated.bim ${out_file}_unrelated.fam

#calculate pcs using flashPCA (need to check what command to use to get 100 pcs; see https://github.com/gabraham/flashpca)
time ./flashpca --bfile ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned -d 100 \
--outload ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_loadings.txt \
--outmeansd ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_meansd.txt 

#write plink files for the whole sample (including related individuals), with only the snps used for pc calculation
time plink \
--bfile ${gen_files} \
--extract ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned.bim \
--make-bed \
--out ${out_file}_maf_0_01_no_long_range_ld_pruned

# project all participants onto the pcs calculated above
time ./flashpca --bfile ${out_file}_maf_0_01_no_long_range_ld_pruned --project --inmeansd ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_meansd.txt --outproj ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_projections.txt --inload ${out_file}_unrelated_maf_0_01_no_long_range_ld_pruned_loadings.txt -v

