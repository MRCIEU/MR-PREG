#!/bin/bash










## Load config.sh 

source config.sh 

## setwd
cd ${pcs}	

module load apps/plink/2.3

plink1='module load apps/plink/1.90'
plink2='module load apps/plink/2.3'

gen_files=${bib_dg}
out_files=${flash}
set_file=price_2008_long_range_ld_snps.txt 


white_mothers=included/white.mums.txt
pakistani_mothers=included/pakis.mums.txt
white_children=included/white.child.txt
pakistani_children=included/pakis.child.txt














## Try loop from old script
#loop over all subgroups set the corresponding outfile and ids to keep
for subgroup in white_mothers pakistani_mothers white_children pakistani_children
#white_mums pakis_mums white_child pakis_child
do 
echo ${subgroup}
out_file=${out_files}${subgroup}
ids=included/"$(cut -d'_' -f1 <<<"$subgroup")"."$(cut -d'_' -f2 <<<"$subgroup")".txt
echo ${ids}
echo ${out_file}

 

 

module load apps/plink/1.90

time plink --bfile ${gen_files} \
 --keep ${ids} \
 --make-set ${set_file} \
 --write-set \
 --threads 10 \
 --out ${out_file}_long_range_ld
 
 
    

module load apps/plink/2.3
 
#remove related individuals
time plink2 \
--bfile ${gen_files} \
--autosome \
--threads 10 \
--keep ${ids} \
--king-cutoff 0.04419417 \
--geno 0.03 \
--make-bed \
--out ${out_file}_unrelated_geno_0_03

#prune snps (flashPCA defaults (https://github.com/gabraham/flashpca), having removed relateds, and remove long range ld
time plink2 \
--bfile ${out_file}_unrelated_geno_0_03 \
--autosome \
--maf 0.01 \
--indep-pairwise 1000 50 0.05 \
--threads 10 \
--geno 0.03 \
--exclude ${out_file}_long_range_ld.set \
--out ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned

module load apps/plink/1.90

#write plink binary files to calculate pcs from
eval $plink1
plink \
--bfile ${out_file}_unrelated_geno_0_03 \
--make-bed \
--extract ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned.prune.in \
--out ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned

#intermediate files can be deleted
rm ${out_file}_unrelated_geno_0_03.bed ${out_file}_unrelated_geno_0_03.bim ${out_file}_unrelated_geno_0_03.fam

#calculate pcs using flashPCA (need to check what command to use to get 100 pcs; see https://github.com/gabraham/flashpca)
time ./flashpca --bfile ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned -d 100 \
--outload ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_loadings.txt \
--outmeansd ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_meansd.txt 

module load apps/plink/2.3

#write plink files for the whole sample (including related individuals), with only the snps used for pc calculation
 
time plink2 \
--bfile ${gen_files} \
--keep ${ids} \
--extract ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned.bim \
--make-bed \
--out ${out_file}_maf_0_01_geno_0_03_no_long_range_ld_pruned

# project all participants onto the pcs calculated above
time ./flashpca --bfile ${out_file}_maf_0_01_geno_0_03_no_long_range_ld_pruned --project --inmeansd ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_meansd.txt --outproj ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_projections.txt --inload ${out_file}_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_loadings.txt -v

done
