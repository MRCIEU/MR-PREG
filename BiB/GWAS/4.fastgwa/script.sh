## GRAMMAR-Gamma is an approximate method that gives O(MN) time of SNP testing instead of O(MN2) 
##(M is the number of SNPs, N - that of individuals)

--fastGWA-mlm --data --est-vg --h2-limit 1.6 --save-fastGWA-mlm-residual --seed 123


#run in the same sample
--fastGWA-lr linear regression
--fastGWA-mlm linear regression
#Generate histogram of pvals to see if the distribution is like Tom’s (normal) or uniform (as expected) – check with GWAS inspector package 
