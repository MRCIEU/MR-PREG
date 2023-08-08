############################################################################
#                                                                          #
#              This script runs MR for each individual study               #
#                                                                          #
############################################################################

############################################################################
#                                Set - UP                                  #
############################################################################

##
## Set working directory 

# from relative path file

##
## Clear the work environment

rm(list = ls())


##
## Setting repository (UoB)

options(repos = c(CRAN ="http://www.stats.bris.ac.uk/R/"))


##
## Setting digits

options(digits = 10)


##
## Library directory

.libPaths()


##
## Updating packages 

#update.packages(ask = "FALSE")


##
## Install packages

#install.packages(c("data.table", "purrr", "devtools", "xlsx", "dplyr"))


##
## Load libraries

library(data.table)
library(dplyr)
library(purrr)
library(haven)
library(TwoSampleMR)
library(MRInstruments)


############################################################################
#                         Load SNP-trait data                              #
############################################################################

##
## Import study specific outcome information

load("./FINAL_METANALYSIS_2021/data/snp-trait.RData")


##
## SNP-BMI from GIANT

giant_expdat <- format_data(giant_dat,
                            type = "exposure",
                            phenotype_col = "exposure",
                            snp_col = "SNP",
                            beta_col = "beta_giant",
                            se_col = "se_giant",
                            eaf_col = "eaf_giant",
                            effect_allele_col = "effect_allele_giant",
                            other_allele_col = "other_allele_giant",
                            pval_col = "pval_giant",
                            samplesize_col = "samplesize_giant"
)


##
## SNP-outcome - study-specific and pooled (FE/RE) + offspring genotype unadjusted / adjusted

all_outdat_tmp <- rename(giant_dat, effect_allele = effect_allele_giant, other_allele = other_allele_giant, eaf = eaf_giant) %>%
  select(., SNP, effect_allele, other_allele, eaf) %>%
  merge(., all_outdat, by.x = "SNP", by.y = "snp") %>%
  filter(., ! (is.na(beta_outcome) | is.na(se_outcome) )) %>%
  mutate(id = paste0(study, "_x_", outcome, "_x_", model)) %>%
  split(., list(.$study, .$model)) 

all_outdat2 <- keep(all_outdat_tmp, ~nrow(.)>0) %>%
  map(., ~format_data(.,
                      type = "outcome",
                      phenotype_col = "outcome",
                      snp_col = "SNP",
                      beta_col = "beta_outcome",
                      se_col = "se_outcome",
                      eaf_col = "eaf",
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      pval_col = "pval_outcome",
                      samplesize_col = "n",
                      id_col = "id"
  ))


############################################################################
#  Run MR and sensitivity analyses (study-specific and pooled estimates)   #
############################################################################

##
## Harmonise SNP-exposure and SNP-outcome data

mr.all_dat <- map(all_outdat2, ~harmonise_data(giant_expdat, ., action = 2)) %>%
  bind_rows %>%
  mutate(
    snp32 = ifelse(SNP %in% snps32$snp, T, F),
    snp97 = ifelse(SNP %in% snps97$snp, T, F)
  )

filter(dat, mr_keep.outcome == F) # 4 observations will be excluded (se = NA)


##
## Extract study, outcome, model for each id

idout <- stringr::str_split(unique(mr.all_dat$id.outcome), "_x_", simplify = TRUE) %>% 
  data.table %>% 
  mutate(id.outcome = unique(mr.all_dat$id.outcome)) %>%
  rename(study = V1, outcome = V2, model = V3)

mr.all_dat <- merge(mr.all_dat, idout[, -"outcome"], by = "id.outcome") 


##
## MR methods

mr.all_res <- mr(mr.all_dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")) %>%
  merge(., idout[, -"outcome"], by = "id.outcome") 

write.csv(mr.all_res, file = "./FINAL_METANALYSIS_2021/data/mr_res.csv")


##
## Sensitivity analyses

fix_dat <- filter(mr.all_dat, study == "Pooled (FE)" & model == "unadjusted")

# Between SNP heterogeneity (ie Cochrane's Q)

het <- mr_heterogeneity(fix_dat, method_list="mr_ivw")

# Directional pleiotropy (ie MR-Egger intercept)

int <- mr_pleiotropy_test(fix_dat)

# Single SNP Wald ratios

single <- mr_singlesnp(fix_dat) %>%
  mutate(SNP = ifelse(SNP=="All - Inverse variance weighted", "All - IVW", SNP))

p_single <- mr_forest_plot(single)


# Leave-one-out

loo <- mr_leaveoneout(fix_dat)

p_loo <- mr_leaveoneout_plot(loo)


##
## save the workspace

rm(list=setdiff(ls(), c("mr.all_dat", "mr.all_res", "het", "int", "loo", "p_loo", "single", "p_single")))

save.image("./FINAL_METANALYSIS_2021/data/mr_res.RData")