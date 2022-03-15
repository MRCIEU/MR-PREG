#remotes::install_github("danjlawson/pcapred.largedata")
#remotes::install_github("danjlawson/pcapred")

library(pcapred)
library(pcapred.largedata)


ref=readreference(ukb_pcs_200())

dat=readbed("/mnt/storage/private/mrcieu/data/bib/genetic/variants/arrays/imputed/uk10k_1000genomes/released/2017-01-27/data/pcs/hardcalled/mergedhc") # Read "your" data
dat=mergeref(dat,ref=ref)     # Merge with the reference (using the included standard reference of 18 UK Biobank Pcs by default)
dat$mode<-"flashpca"
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
writepred("../../data/pcs/output/pcapred/projected.eigenvals_200",dat,pred) # Write output in plink --covar format