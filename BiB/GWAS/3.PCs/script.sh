#Adapted to bib data by Marwa Al Arab, based on script from Tom Bond.
#A script to i)generate pcs for 4 stratas in bib data: White mothers, Pakistanis mothers, White Children & Pakistanis Children
#!/bin/bash

#SBATCH --job-name=pcs_job
#SBATCH --partition=mrcieu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --array=1-4


#calculate pcs for unrelated individuals, then project all individuals including related individuals onto these pcs
#change directory to scripts
cd $bibimputed/scripts/pcs/

#define variables 
plink1='module load apps/plink/1.90'
plink2='module load apps/plink/2.00'
eval $bibimputed
eval $bibdirect
#directly genotyped data 
gen_files=$bibdirect/data
#hardcalled_files = ../../data/pcs/hardcalled/chr${SLURM_ARRAY_TASK_ID}
out_files=$bibimputed/scripts/pcs/data/pcs/output/
#out_path=
set_file=price_2008_long_range_ld_snps.txt 
white_mums_ids=included/white.mums.txt
pakis_mums_ids=included/pakis.mums.txt
white_child_ids=included/white.child.txt
pakis_child_ids=included/pakis.child.txt

case ${SLURM_ARRAY_TASK_ID} in 
1)
ids=white_mums_ids
;;
2)
ids=pakis_mums_ids
;;
3)
ids=white_child_ids
;;
4)
ids=pakis_child_ids
;;
*)
echo "Task Array id not matching case!!"
esac


out_file=${out_files}ids




#This command is to generate hardcalled files to use for merging data with UKB pcs with pcapred
#to write a set of plink bed/bim/fam files, which contain “hard called” genotype at imputed snps (i.e. integer values representing the most likely imputed genotype).
#module load apps/plink/2.00
#time plink --vcf $gen_file dosage=GP --mach-r2-filter 0.82 --hard-call-threshold 0.1 --import-dosage-certainty 0.9 --make-bed --out $hardcalled_file
 
eval $plink1
#module load apps/plink/1.90
#make .set file with a list of variants in long range ld (from Price, Alkes L., et al. "Long-range LD can confound genome scans in admixed populations." The American Journal of Human Genetics 83.1 (2008): 132-135.)
 time plink --bfile ${gen_files} \
 --keep ${ids} \
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
--keep ${ids} \
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

