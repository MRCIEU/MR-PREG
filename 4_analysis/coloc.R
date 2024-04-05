############################################################################
#                                Set - UP                                  #
############################################################################

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

#install.packages(c("data.table", "purrr", "devtools", "xlsx", "readxl", "dplyr", "ggplot2", "remotes"))
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments")
#devtools::install_github("NightingaleHealth/ggforestplot")
#remotes::install_github("mattlee821/EpiCircos")


##
## Load libraries

library(data.table)
library(purrr)
library(TwoSampleMR)
library(MRInstruments)
#library(xlsx)
library(readxl)
library(dplyr)


##
## List available GWASs in MR-Base catalogue

ao <- available_outcomes()


##
## Exposure info file

exp.info <- read_excel("data/exp_info.xlsx", sheet = "iv_exp")


##
## NMR metabolites info file

nmr.info <- read_excel("data/biomarker_annotations_updated (1).xlsx")
  

############################################################################
#                       Prepare and harmonise data                         #
############################################################################

##
## Exposure data

expdat <- read.table("data/expdat.txt")


##
## Outcome data

outdat <- read.table("data/nmr/nmr_female_sumdat_ivs.txt") %>%
    			rename_all(., .funs = list(~ tolower(.))) %>%
    			rename(effect_allele = allele1, other_allele = allele0, eaf = a1freq, pval = p_bolt_lmm_inf) %>%
				filter(info >= 0.8) %>%
    			format_data(.,
    							type = "outcome",
    							phenotype_col = "id",
    							snp_col = "snp",
    							pos_col = "bp"
    						)
 
length(unique(outdat$outcome)) # Observed number of outcomes (expected: 249 metabolic measures)


##
## Harmonise exposure and outcome data

dat <- harmonise_data(expdat, outdat)


############################################################################
#                                  Run MR                                  #
############################################################################

##
## Run MR

res.raw <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))			

 
##
## Format MR results

res <-  res.raw %>%
				merge(., exp.info, by.x = "id.exposure", by.y = "id") %>%
				mutate(outcome = gsub("nmr_dat_female_|_int_imputed.txt.gz_TMPDAT.txt", "", outcome)) %>%
				merge(., nmr.info, by.x = "outcome", by.y = "Name in TSV",  all.x = T) %>%
				rename(Derived = "Derived measure", metab_name = "Biomarker name", priority = analysis) %>%
				mutate(
						Category = ifelse(Group %in% c("Lipoprotein subclasses", "Relative lipoprotein lipid concentrations"), as.character(Subgroup), as.character(Group)),
						Derived = as.character(Derived),
						method_abbrev = case_when(method == "Inverse variance weighted" ~ "ivw",
												  method == "MR Egger" ~ "mregger",
												  method == "Weighted median" ~ "wmedian"
					   ),
					   method = factor(method, levels=c("MR Egger", "Weighted median", "Inverse variance weighted")),
					   type_analyses = "MR",
					   assay = "nmr"
					  )
			

##
## Save MR results

mr.res <- res %>%
        	arrange(trait, Category, metab_name, method) %>%
        	select(trait, abbrev, id.exposure, id.outcome, outcome, metab_name, Category, Derived, assay, nsnp, type_analyses, priority, method, b, se, pval) 
        	
writexl::write_xlsx(list(mr_res = mr.res), path = "./results/mr_results.xlsx")

				
q("no")					