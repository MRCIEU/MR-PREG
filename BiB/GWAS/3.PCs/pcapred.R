#Marwa Al Arab
#generate UKB pcs after merging the local data (bib) with reference data (see https://github.com/danjlawson/pcapred for documentation)
#remotes::install_github("danjlawson/pcapred.largedata")
#remotes::install_github("danjlawson/pcapred")

library(pcapred)
library(pcapred.largedata)
#get the location of the data from linux environment variables
imputed_data<-(getenv "bibimputed")

ref=readreference(ukb_pcs_200())

dat=readbed(paste0(imputed_data,"data/pcs/hardcalled/mergedhc")) # Read "your" data
dat=mergeref(dat,ref=ref)     # Merge with the reference (using the included standard reference of 18 UK Biobank Pcs by default)
dat$mode<-"flashpca"
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
writepred(paste0(imputed_data,"data/pcs/output/pcapred/projected.eigenvals_200"),dat,pred) # Write output in plink --covar format