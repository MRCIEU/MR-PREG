############################################################################
#                                                                          #
#  This script creates plots sensitivity MV (pre/during weight) analyses   #
#                                                                          #
############################################################################

############################################################################
#                                Set - UP                                  #
############################################################################

##
## Set working directory 



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
#library(forestplot)
library(ggplot2)
library(ggforestplot)
#library(haven)
#library(TwoSampleMR)
#library(MRInstruments)
library(tidyverse)

############################################################################
#                             Format data                                  #
############################################################################

##
## Load data

# MV 
load("./FINAL_METANALYSIS_2021/data/mv_res_wgt.RData")


##
## Combine datasets

out_info2 <-  filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB)) 

mv.all_res_wgt <- mutate(mv.all_res, analyses = "MV") %>%
  mutate(outname2 = ifelse(outcome == "gestational_diabetes", "Gestational diabetes", outname))

mv.all_res_wgt <- mv.all_res_wgt %>% drop_na(ind)

mv.all_res_wgt <- subset(mv.all_res_wgt, study == "Pooled (FE)")

## Directory for saving outputs

out_dir <- "FINAL_METANALYSIS_2021/output/"

## Study-specific and pooled (fixed- and random-effects) IVW results 
## Function to select data subset

sel_dat <- function(dat, studies, models, typevar) {
  
  # Filter relevant data  
  df <- filter(dat, 
               ind %in% studies
               &
                 model %in% models 
               & 
                 type %in% typevar
  ) %>%
    mutate(      group = factor(group, levels = c("Pregnancy", "Delivery", "Birth", "Postnatal"))
      #,
      # study = factor(study, levels=rev(levels("ALSPAC", "BIB", "DNBC", "EFSOCH", "FinnGen", "GEN3G", "GENR", "GOYA", "HAPO", "INMA", 
      #                                        "MOBA", "NFBC1966", "NFBC1986",  "UK_BIOBANK", "Pooled (FE)",  "Pooled (RE)")))
    ) %>%
    arrange(order) 
}           


##
## Forest plot function

fplot_dat <- function(dat, strata, or) {
  
  dat  %>%
    ggforestplot::forestplot(
      df = .,
      name = outname2,
      estimate = b,
      pvalue = pval,
      logodds = or,
      psignif = 0.05,
      xlab = "OR (95% CI)\nper 1-SD increment in BMI",
      colour = strata
    ) +
    theme(
      legend.text = element_text(size = 15),
      legend.title = element_blank(),
      axis.text.x = element_text(size=15),
      axis.title.x = element_text(size=15),
      axis.text.y = element_text(size=15),
      strip.text = element_text(size = 17)
    ) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  
}            

# Binary traits 

##pooled effects
## Pooled (fixed- and random-effects) IVW results 


df3a <- sel_dat(mv.all_res_wgt, c("During pregnancy", "Pre-pregnancy"), "adjusted", "binary")

fplot_dat(df3a, df3a$ind, or = T) 
.
ggsave(paste0(out_dir, "forest_mv-pooleff_bin_suppfig_wgt.jpg"),  width = 15, height = 8)


# Continuous traits 

df3b <- sel_dat(mv.all_res_wgt, c("During pregnancy", "Pre-pregnancy"), "adjusted", "continuous")

fplot_dat(df3b, df3b$ind, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mv-pooleff_cont_suppfig_wgt.jpg"),  width = 15, height = 8) 





