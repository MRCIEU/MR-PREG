############################################################################
#                                                                          #
#       This script creates sensitivity plot comparing PNC models          #
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

library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MRInstruments)
library(purrr)
library(openxlsx)
library(meta)
library(metafor)
library(cowplot)
library(gridGraphics)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(readxl)
library(data.table)
library(ggforestplot)


############################################################################
#                             Format data                                  #
############################################################################

##
## Load data

MV_main1_c <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "Continuous Adjusted")
MV_main1_c$models <-"A Main analysis"
MV_main1_b <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "Binary Adjusted")
MV_main1_b$models <-"A Main analysis"

MV_main2_c <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "2 Maternal BMI CONTINUOUS")
MV_main2_c$models <-"B Main analysis, restricted to PNC sample"
MV_main2_b <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "2 Maternal BMI BINARY")
MV_main2_b$models <-"B Main analysis, restricted to PNC sample"

MV_main3_c <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "3 Maternal BMI CONTINUOUS")
MV_main3_c$models <-"C Main analysis, adjusted for paternal BMI"
MV_main3_b <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "3 Maternal BMI BINARY")
MV_main3_b$models <-"C Main analysis, adjusted for paternal BMI"

MV_main4_c <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "2 Paternal BMI CONTINUOUS")
MV_main4_c$models <-"D Paternal BMI, adjusted for maternal BMI"
MV_main4_b <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis.xls" , sheet = "2 Paternal BMI BINARY")
MV_main4_b$models <-"D Paternal BMI, adjusted for maternal BMI"

#put mv_main together

MV_main_b <- bind_rows(MV_main1_b, MV_main2_b, MV_main3_b, MV_main4_b)
MV_main_c <- bind_rows(MV_main1_c, MV_main2_c, MV_main3_c, MV_main4_c)


#---------------------------------------------------------------------#
#               meta-analyse mat and PNC                              #
#---------------------------------------------------------------------#
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")
out_info2 <- filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB))

mv.all_res <- mutate(mv.all_res, analyses = "MV") 
#BINARY
#merge with Carolinas outcomes
MV_main_b$outcome[MV_main_b$outcome=="sga"] <- 'sga_all'
MV_main_b$outcome[MV_main_b$outcome=="lga"] <- 'lga_all'

names(MV_main_b)[names(MV_main_b) == 'outcome'] <- 'phase2_varname'

res_MV_b <- merge(MV_main_b, out_info, by = "phase2_varname") 

names(res_MV_b)[names(res_MV_b) == 'logOR'] <- 'beta'

as.factor(res_MV_b$group)

res_MV_b <- arrange(res_MV_b, models)
as.factor(res_MV_b$models)

res_MV_b$models <- factor(res_MV_b$models, levels = c("D Paternal BMI, adjusted for maternal BMI", "C Main analysis, adjusted for paternal BMI", "B Main analysis, restricted to PNC sample", "A Main analysis"))
res_MV_b <- arrange(res_MV_b, order, group, models)

#CONTINUOUS
#merge with Carolinas outcomes
MV_main_c$outcome[MV_main_c$outcome=="ponderal_zscore"] <- 'ponderal_index'
MV_main_c$outcome[MV_main_c$outcome=="birthlength_zscore"] <- 'birthlength'
MV_main_c$outcome[MV_main_c$outcome=="zbw_subsamp"] <- 'birthweight'
MV_main_c$outcome[MV_main_c$outcome=="gestage_zscore"] <- 'gestational_age'

names(MV_main_c)[names(MV_main_c) == 'outcome'] <- 'phase1_varname'

res_MV_c <- merge(MV_main_c, out_info, by = "phase1_varname") 


as.factor(res_MV_c$group)

res_MV_c <- arrange(res_MV_c, models)
as.factor(res_MV_c$models)

res_MV_c$models <- factor(res_MV_c$models, levels = c("D Paternal BMI, adjusted for maternal BMI", "C Main analysis, adjusted for paternal BMI", "B Main analysis, restricted to PNC sample", "A Main analysis"))
res_MV_c <- arrange(res_MV_c, order, group, models)


## Directory for saving outputs

out_dir <- "FINAL_METANALYSIS_2021/output/"

## Study-specific and pooled (fixed- and random-effects) IVW results 
## Function to select data subset

sel_dat <- function(dat, typevar) {
  
  # Filter relevant data  
  df <- filter(dat, 
                 type %in% typevar
  ) %>%
    mutate(      group = factor(group, levels = c("Pregnancy", "Delivery", "Birth", "Postnatal")),
      models = factor(models, levels=c("D Paternal BMI, adjusted for maternal BMI", "C Main analysis, adjusted for paternal BMI",  "B Main analysis, restricted to PNC sample", "A Main analysis"))
    ) %>%
    arrange(order, models) 
}           


##
## Forest plot function
# Set class to factor to set order of display.

fplot_dat <- function(dat, or) {
  
  dat  %>%    ggforestplot::forestplot(
      df = .,
      name = outname,
      estimate = beta,
      pvalue = pval,
      logodds = or,
      psignif = 0.05,
      xlab = "OR (95% CI)\nper 1-SD increment in BMI",
      colour = models
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

# Binary traits - part 1

##pooled effects - by all models
df3a <- sel_dat(res_MV_b, "binary")

fplot_dat(df3a, or = T) 

ggsave(paste0(out_dir, "forest_pnc_bin_allmod_suppfig_ALSPAC.jpg"),  width = 15, height = 8)


df3b <- sel_dat(res_MV_c , "continuous")

fplot_dat(df3b, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_pnc_cont_allmod_suppfig_ALSPAC.jpg"),  width = 15, height = 8) 




