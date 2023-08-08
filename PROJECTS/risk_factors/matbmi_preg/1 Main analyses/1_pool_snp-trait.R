############################################################################
#                                                                          #
#         This script pools SNP-trait association data across studies      #
#                                                                          #
############################################################################

############################################################################
#                                Set - UP                                  #
############################################################################

##
## Set working directory 

#from relative path file


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
library(meta)
#library(metafor)
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)


############################################################################
#              Reference file for BMI-increasing alleles                   #
############################################################################

##
## List of IVs (32 and 97 SNP sets)

snps32 <- read_dta("./META_ANALYSIS/descriptives/32_riskscore.dta") %>% 
  mutate(set = 32) %>% 
  rename_all(., ~sub("_32", "", .)) 

snps97 <- read_dta("./META_ANALYSIS/descriptives/97_riskscore.dta") %>% 
  mutate(set = 97) %>% 
  rename_all(., ~sub("_97", "", .))

snps <- bind_rows(snps32, snps97)

snps_unique <- c(snps32$snp, snps97$snp) %>% unique


##
## Import SNP-BMI data from GIANT European females

giant_dat <- extract_outcome_data(snps = snps_unique, outcomes = 'ieu-a-974')  %>%
  rename_all(., ~sub(".outcome", "_giant", .)) %>%
  mutate(
    eaf_exposure = ifelse(beta_giant >= 0, eaf_giant, (1-eaf_giant) ),
    ea = ifelse(beta_giant >= 0, effect_allele_giant, other_allele_giant),
    other_allele_giant = ifelse(beta_giant >= 0, other_allele_giant, effect_allele_giant),
    effect_allele_giant = ea,
    beta_giant = abs(beta_giant),
    exposure = "BMI (GIANT females)"
  ) %>%
  select(-ea, -outcome, -outcome.deprecated)

ea_97snps <- merge(giant_dat, snps97, by.x = "SNP", by.y = "snp")
table(ea_97snps$ea, ea_97snps$effect_allele_giant) # All EA match between GIANT and MR-PREG studies

ea_32snps <- merge(giant_dat, snps32, by.x = "SNP", by.y = "snp")
table(ea_32snps$ea, ea_32snps$effect_allele_giant)  # All EA match between GIANT and MR-PREG studies


############################################################################
#                      Import auxiliary files                              #
############################################################################

##
## Load study info file (studies and set of 32 x 97 SNPs)

studies <- read.csv("./FINAL_METANALYSIS_2021/data/list_of_studies.csv")

studies_p1 <- filter(studies, phase1 == 1 & keep == "yes") %>% pull(study_var)


##
## Load outcome info file 

out_info <- read.csv("./FINAL_METANALYSIS_2021/data/list_of_varnames.csv") %>%
 # rename(outname = "?..Outcome") %>%
  rename(outname = Outcome) %>%
  filter(! is.na(outname)) %>%
  filter(primary == "yes")


############################################################################
#                   Import and format files - phase 1                      #
############################################################################

##
## Extract summary data for studies from phase 1

expvar <- c("study", "outcome", "snp", "n", "beta_exposure", "se_exposure", "r2_exposure", "fstat_exposure")

outvar <- c("study", "outcome", "snp", "n", "beta_outcome", "se_outcome", "pval_outcome")

extr_dat <- function(study, vars, model = NULL) {
  
  sumdat <- read_dta(paste0("./META_ANALYSIS/single_study_results/", study, "_snp_results", model, ".dta")) %>%
    rename_all(., ~tolower(.)) 
  
  if("p_outcome" %in% names(sumdat)) {
    
    sumdat <- rename(sumdat, pval_outcome = p_outcome)
    
  }
  
  dat <- mutate(sumdat, study = toupper(study)) %>%
    select(., all_of(vars))
} 


##
## SNP-exposure estimates (not available for UKBB)
expdat <- map(studies_p1[studies_p1 != "UK_BIOBANK"], ~extr_dat(., expvar)) %>% 
  bind_rows %>%
  mutate(exposure = "Maternal BMI")


##
## SNP-outcome estimates (not adjusted for offspring genotype)
outdat_p1 <- map(studies_p1, ~extr_dat(., outvar)) %>% 
  bind_rows %>%
  mutate(outcome = ifelse(outcome=="macro", "macrosomia", outcome),
         outcome = ifelse(outcome=="membrane", "membrane_rupture", outcome)
  ) %>%
  filter(! study %in% c("ALSPAC", "UK_BIOBANK")) %>%
  mutate(model = "unadjusted")


##
## SNP-outcome estimates (adjusted for offspring genotype) - only available for 6 studies in phase1
outdat_p1_adj <- filter(studies, offspring_adjusted == "yes" & phase1 == 1) %>%
  pull(study_var) %>%
  map(., ~extr_dat(., outvar, model="_adj")) %>% 
  bind_rows %>% 
  mutate(outcome = ifelse(outcome=="macro", "macrosomia", outcome),
         outcome = ifelse(outcome=="membrane", "membrane_rupture", outcome)
  ) %>%
  mutate(model = "adjusted")


############################################################################
#                   Import and format files - phase 2                      #
############################################################################

##
## Studies from phase 2

studies_p2 <- filter(studies, phase2 == 1 & keep == "yes") %>% pull(study_var)


##
## Extract summary data for studies from phase 2

p94_dir <- "//rdsfcifs.acrc.bris.ac.uk/MRC-IEU-research/projects/ieu2/p6/094/working/results/"

extr_dat2 <- function(source, a, b) {
  
  # Harmonise outcome name within phase 2 and with phase 1 studies
  if(source == "FinnGen") {
    
    out_info <- mutate(out_info, phase2_varname2 = phase2_varname_Finn)
    
  } else if (source == "UK_BIOBANK") {
    
    out_info <- mutate(out_info, phase2_varname2 = phase2_varname_UKBB)
    
  } else {
    
    out_info <- mutate(out_info, phase2_varname2 = phase2_varname)
    
  }
  
  outvar_info <- filter(out_info, ! (phase1_varname == "" |  phase2_varname2 == "")) %>%
    select(outname, phase1_varname, phase2_varname2) 
  
  outvar <- outvar_info$phase2_varname2
  
  # Import data, keep BMI SNPs and format var names
  sumdat <- map(outvar, ~fread(paste0( p94_dir, "GWAS/", a, ., b))) %>%
    map(., ~rename_all(., ~tolower(.))) %>%
    map(., ~filter(., snp %in% snps97$snp)) %>%
    map2(., outvar, ~mutate(.x, study = source, outcome = .y)) %>%
    bind_rows
  
  if("p" %in% names(sumdat)) {
    sumdat <- rename(sumdat, pval = p)
  }
  
  sumdat <- sumdat %>%
    rename(beta_outcome = beta, se_outcome = se, pval_outcome = pval) 
  
  # Harmonise data so that effect allele = BMI increasing allele
  sumdat2 <- merge(sumdat, giant_dat, by.x = "snp", by.y = "SNP") %>%
    mutate(beta_outcome = ifelse(effect_allele == effect_allele_giant, beta_outcome, -beta_outcome)) %>%
    merge(.,  outvar_info, by.x = "outcome", by.y = "phase2_varname2") %>%
    select(-outcome) %>%
    mutate(outcome = phase1_varname)
  
} 

memory.limit(size=250000)

##
## SNP-outcome estimates (not adjusted for offspring genotype)

# ALSPAC
alspac_dat <- extr_dat2(source = "ALSPAC", a = "ALSPAC/mothers/ALSPACmums.", b = ".20210716.txt.gz") %>%
  select("study", "outcome", "snp", "beta_outcome", "se_outcome", "pval_outcome", "samplesize", "ncase", "ncontrol") %>%
  rename(n = samplesize)

# MOBA 
moba_dat <- extr_dat2(source = "MOBA", a = "MOBA-17k/mothers/k17_mum_comb_", b = ".txt") %>%
  select("study", "outcome", "snp", "beta_outcome", "se_outcome", "pval_outcome", "n_total", "n_cases", "n_controls") %>%
  rename(n = n_total, ncase = n_cases, ncontrol = n_controls)

# FinnGen 
finn_dat <- extr_dat2(source = "FinnGen", a = "FINNGEN/mothers/FINNGEN-R5.", b = ".20210812.txt.gz") %>%
  select("study", "outcome", "snp", "beta_outcome", "se_outcome", "pval_outcome", "samplesize", "ncase", "ncontrol")  %>%
  rename(n = samplesize)

# UKBB - NEED TO ADD SAMPLE SIZE TO UKBB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ukb_dat <- extr_dat2(source = "UK_BIOBANK", a = "UKB/mothers/", b = ".plink.gz") %>%
  select("study", "outcome", "snp", "beta_outcome", "se_outcome", "pval_outcome") %>%
  distinct(.keep_all = T)

# Combined
outdat_p2 <- bind_rows(alspac_dat, moba_dat, finn_dat, ukb_dat) %>%
  mutate(model = "unadjusted")


##
## SNP-outcome estimates (adjusted for offspring genotype) - only available for ALSPAC and MOBA

#st <- c("ALSPAC", "MOBA")
st <- c("MOBA")

#st_f <- c ("alspac", "moba17k")
st_f <- c("moba17k")

t <- filter(out_info, ! (phase2_varname=="" | group == "Pregnancy") ) %>% pull(phase2_varname) ## Pregnancy variable are not eligible for offspring adjusted analyses

outdat_p2_adj <- map2(st, st_f, ~fread(paste0(p94_dir, "adjusted-ipd/", .x, "/", .y, "_adj-ipd_sumdat.txt"))) %>%
  map2(st, ., ~mutate(.y, study = .x, model = "adjusted")) %>%
  bind_rows %>%
  filter(
    snp %in% snps97$snp
    &
      effect == "maternal"
    &
      type == "adjusted"
    &
      trait %in% t
  ) %>%
  merge(.,  out_info, by.x = "trait", by.y = "phase2_varname") %>%
  rename(outcome = phase1_varname, beta_outcome = beta, se_outcome = se, pval_outcome = pval) %>%
  select("study", "outcome", "snp", "beta_outcome", "se_outcome", "pval_outcome", "model") 


############################################################################
#                     Metanalyse SNP-exposure associations                 #
############################################################################

##
## Subset data

expdat2 <- group_by(expdat, study, snp) %>% # Group by study
  top_n(1, n) %>% # Keep unique SNP-BMI associations (choose outcome with maximum N) 
  select(-outcome) %>%
  distinct(study, snp, .keep_all = T) 
#filter(expdat, snp %in% snps97$snp) %>% # Keep SNPs in the 97 SNP set         

##
## Function to perform metanalysis of SNP-trait associations across studies

run_meta_bmi <- function(rsid) {
  
  # Keep SNP-BMI data for one SNP at a time
  input <- filter(expdat2, snp == rsid) 
  print(input)
  
  # Meta-analysis across studies
  ma <- metagen(TE = beta_exposure, 
                seTE = se_exposure, 
                data = input,
                studlab = study, 
                backtransf = F,
                hakn = F, 
                method.tau="DL", 
                comb.fixed = T, 
                comb.random = T
  )
  print(ma)
  
  # Extract values from fixed-effect metanalysis
  fix_ma <- data.table(
    beta_exposure = ma$TE.fixed,
    se_exposure   = ma$seTE.fixed,
    pval_exposure = ma$pval.fixed,
    nstudies      = ma$k,
    q             = ma$Q,
    q_pval        = ma$pval.Q,
    i2            = ma$I2,
    exposure      = "Maternal BMI",
    SNP           = rsid,
    study         = "Pooled (FE)"
  )
  
  return(fix_ma)
}


##
## Metanalyse SNP-BMI estimates across studies - Need to add MOBA SNP-BMI estimates!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pooled_expdat <- map(snps_unique, ~run_meta_bmi(.)) %>% 
  bind_rows %>%
  mutate(set = case_when(
    SNP %in% intersect(snps32$snp, snps97$snp) ~ "Both sets",
    SNP %in% snps97$snp ~ "97 SNPs set",
    SNP %in% snps32$snp ~ "32 SNPs set"
    
  )) 
table(pooled_expdat$set)
# 32 SNPs set 97 SNPs set   Both sets 
# 20          85          12 


##
## Scatter plot GIANT females x Maternal BMI

all_expdat <- merge(pooled_expdat, giant_dat, by = "SNP") %>%
  mutate(ll_ci_exposure = beta_exposure - 1.96 * (se_exposure), 
         ul_ci_exposure = beta_exposure + 1.96 * (se_exposure),
         ll_ci_giant    = beta_giant - 1.96 * (se_giant), 
         ul_ci_giant    = beta_giant + 1.96 * (se_giant),
         r2_exposure    = round(2 * beta_exposure ^ 2 * eaf_exposure * (1-eaf_exposure), 3),
         r2_giant       = round(2 * beta_giant ^ 2 * eaf_exposure * (1-eaf_giant), 3),
         f_exposure     = round( ((beta_exposure/se_exposure)^2), 1),
         f_giant        = round( ((beta_giant/se_giant)^2), 1)
  )
# Total R2
filter(all_expdat, SNP %in% snps97$snp) %>%
  summarise(R2_exposure = sum(r2_exposure), R2_giant = sum(r2_giant))
#R2_exposure R2_giant
#1       0.016    0.027

# Total R2
filter(all_expdat, SNP %in% snps97$snp) %>%
  summarise(F_exposure = mean(f_exposure), F_giant = mean(f_giant))
#   F_exposure   F_giant
#1     4.2625 36.271875

# Correlation SNP-BMI between GIANT and cohorts
bmi_cor <- filter(all_expdat, SNP %in% snps97$snp)
cor(bmi_cor$beta_exposure, bmi_cor$beta_giant, use = "complete.obs")
# [1] 0.6660248015  

p1 <-  ggplot(all_expdat, aes(x = beta_giant, y = beta_exposure)) + 
  geom_errorbar(aes(ymin = ll_ci_exposure, ymax = ul_ci_exposure), color = "lightgrey", alpha = 1.5, size = 0.3) +
  geom_errorbarh(aes(xmin = ll_ci_giant, xmax = ul_ci_giant), color = "lightgrey", alpha = 1.5, size = 0.3) +
  geom_point(aes(color = set, shape = set), size = 3.5, alpha = 0.6) +
  theme_classic() +
  scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size=16)
  ) +
  scale_x_continuous(name = "BMI (GIANT females)", limits=c(-0.002, 0.1), breaks=seq(0,0.1,0.025)) +
  scale_y_continuous(name = "Maternal BMI", limits=c(-0.11, 0.18), breaks=seq(-0.1,0.15, 0.05)) +
  geom_abline(slope=1, intercept=0, color="grey")

ggsave("./FINAL_METANALYSIS_2021/output/scatter_snp-bmi.pdf", width = 20, height = 20, units = "cm")

write.csv(all_expdat, file = "./FINAL_METANALYSIS_2021/data/snp-bmi.csv")


############################################################################
#             Metanalyse SNP-outcome associations (97 SNPs set only)       #
############################################################################

##
## Combine data from phase 1 and phase 2 studies

outdat <-  bind_rows(outdat_p1, outdat_p1_adj) %>%
  bind_rows(., outdat_p2, outdat_p2_adj)


##
## Function to perform metanalysis of SNP-outcome associations across studies

run_meta <- function(trait, rsid, m) {
  
  # Subset dataset to keep estimates for one trait, one SNP, one model across studies
  input <- filter(outdat, outcome == trait & snp == rsid & model == m)
  print(input)
  
  if(nrow(input) == 0) {
    
    return(NULL)
    
  } else {
    
    # Meta-analysis across studies
    ma <- metagen(TE = beta_outcome, 
                  seTE = se_outcome, 
                  data = input,
                  studlab = study, 
                  backtransf = F,
                  hakn = F, 
                  method.tau="DL", 
                  comb.fixed = T, 
                  comb.random = T
    )
    print(ma)
    
    # Extract values from fixed-effect metanalysis
    fix_ma <- data.table(
      beta_outcome = ma$TE.fixed,
      se_outcome   = ma$seTE.fixed,
      pval_outcome = ma$pval.fixed,
      nstudies     = ma$k,
      q            = ma$Q,
      q_pval       = ma$pval.Q,
      i2           = ma$I2,
      study        = "Pooled (FE)"
    )
    
    # Extract values from random-effect metanalysis
    ran_ma <- data.table(
      beta_outcome = ma$TE.random,
      se_outcome   = ma$seTE.random,
      pval_outcome = ma$pval.random,
      nstudies     = ma$k,
      q            = ma$Q,
      q_pval       = ma$pval.Q,
      i2           = ma$I2,
      study        = "Pooled (RE)"
    )
    
    # Combined results from FE and RE metanalysis
    all_ma <- bind_rows(fix_ma, ran_ma) %>%
      mutate(outcome = trait, snp = rsid, model = m)
    
    return(all_ma)
    
  }
  
}


##
## Arguments for metanalysis function 

args_meta <- cross_df(list(
  trait = unique(outdat$outcome), 
  rsid  = snps97$snp,
  m = unique(outdat$model)
))


##
## Run metanalyses

meta_outdat <- pmap(args_meta, run_meta) %>% bind_rows


##
## Combine all outcome data

all_outdat <- bind_rows(outdat, meta_outdat)


##
## save the workspace

rm(list=setdiff(ls(), c("giant_dat", "out_info", "all_outdat", "snps32", "snps97")))

save.image(file = "./FINAL_METANALYSIS_2021/data/snp-trait.RData")


