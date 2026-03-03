###############################################################################
########### script to combine array and PCs 1-40 to create ####################
########### covariate files to run UKB GWAS using REGENIE #####################
###############################################################################I







### set up ###################################################################
library(dplyr)
rm(list = ls())

scratch <- Sys.getenv("MRPREG_sumdat")
setwd(paste0(scratch, ""))

# now using latest 2018 imputed path -
ukb <- Sys.getenv("UKB_DIR")
ukb_dir <- paste0(ukb, "/released/2018-09-18/data/derived/")

out_dir <- "./data/GWAS/UKB/mothers/"

### read in covariate files ###################################################
sex_array <- read.table(paste0(ukb_dir, "standard_covariates/data.covariates.plink.txt"),
    header = FALSE)
colnames(sex_array) <- c("FID", "IID", "sex", "array")
# check
dim(sex_array); head(sex_array)

PCs <- read.table(paste0(ukb_dir, "principal_components/data.pca1-40.plink.txt"),
    header = FALSE)
colnames(PCs) <- c("FID", "IID", paste0("PC", seq(1:40)))
# check
dim(PCs); head(PCs)

# remove sex column
array <- sex_array[,c(-3)]
# change array to 1 & 2
array <- array %>%
    mutate(array = case_when(array == "UKBB" ~ 0,
                        array == "UKBL" ~ 1))
# check
dim(array); head(array)

### combine ###################################################################
covs <- left_join(array, PCs, by = c("FID", "IID"))
# check
dim(covs); head(covs)

### write out #################################################################
write.table(covs,
    paste0(out_dir, "covs.txt"),
    quote = FALSE, row.names = FALSE)
