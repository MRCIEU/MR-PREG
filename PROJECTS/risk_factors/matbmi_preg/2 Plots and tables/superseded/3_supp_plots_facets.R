############################################################################
#                                                                          #
#  This script creates plots/tables for main and sensitivity MR analyses   #
#                                                                          #
############################################################################

############################################################################
#                                Set - UP                                  #
############################################################################

##
## Set working directory 

setwd("//rdsfcifs.acrc.bris.ac.uk/MRC-IEU-research/projects/ieu2/p6/011/working/data/analysis/")


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


############################################################################
#                             Format data                                  #
############################################################################

##
## Load data

# MV 
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")

# MR
load("./FINAL_METANALYSIS_2021/data/mr_res.RData")


##
## Combine datasets

out_info2 <-  filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB)) 

mv.all_res <- mutate(mv.all_res, analyses = "MV") %>%
  mutate(outname2 = ifelse(outcome == "gestational_diabetes", "Gestational diabetes", outname))

mr.all_res <- mutate(mr.all_res, analyses = "MR") %>%
  merge(., out_info2, by.x = "outcome", by.y = "phase1_varname") %>%
  mutate(outname2 = ifelse(outcome == "gestational_diabetes", "Gestational diabetes", outname))

all_res <- bind_rows(mv.all_res, mr.all_res) 

## Directory for saving outputs

out_dir <- "FINAL_METANALYSIS_2021/output/"

## Study-specific and pooled (fixed- and random-effects) IVW results 
## Function to select data subset

## Study-specific colour scheme (14 studies)

color <- c("aquamarine2", "black", "blue2", "blueviolet", "brown4", "darkcyan", 
           "chartreuse", "coral2", "darkgoldenrod3", "darkmagenta", "darkolivegreen4", "deeppink2", 
           "deepskyblue1", "burlywood3"
)


sel_dat <- function(dat, studies, models, typevar) {
  
  # Filter relevant data  
  df <- filter(dat, 
               study %in% studies
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

# Binary traits - part 1

df1a <- sel_dat(mv.all_res, unique(all_res$study), "adjusted", "binary") %>% filter(group %in% c("Pregnancy", "Delivery"))

fplot_dat(df1a, df1a$study, or = T) + scale_color_manual(values=color)

ggsave(paste0(out_dir, "forest_mv-studies_bin_part1_suppfig4A2.jpg"),  width = 15, height = 10)

# Binary traits - part 1

df1b <- sel_dat(mv.all_res, unique(all_res$study), "adjusted", "binary") %>% filter(group %in% c("Birth", "Postnatal"))

fplot_dat(df1b, df1b$study, or = T) + scale_color_manual(values=color)

ggsave(paste0(out_dir, "forest_mv-studies_bin_part2_suppfig4A2.jpg"),  width = 15, height = 10)

# Continuous traits 

df1c <- sel_dat(mv.all_res, unique(all_res$study), "adjusted", "continuous")

fplot_dat(df1c, df1c$study, or = F) + xlab("MD (95% CI)\nper 1-SD increment in BMI") + scale_color_manual(values=color)

ggsave(paste0(out_dir, "forest_mv-studies_cont_suppfig4B2.jpg"),  width = 15, height = 10)


##pooled effects
## Pooled (fixed- and random-effects) IVW results 

df3a <- sel_dat(mv.all_res, c("Pooled (FE)", "Pooled (RE)"), "adjusted", "binary")

fplot_dat(df3a, df3a$study, or = T) 

ggsave(paste0(out_dir, "forest_mv-pooleff_bin_suppfig5A.jpg"),  width = 15, height = 8)


df3b <- sel_dat(mv.all_res, c("Pooled (FE)", "Pooled (RE)"), "adjusted", "continuous")

fplot_dat(df3b, df3b$study, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mv-pooleff_cont_suppfig5B.jpg"),  width = 15, height = 8) 


##pooled effects - by all models
df3a <- sel_dat(mv.all_res, c("Pooled (FE)"), c("basic", "adjusted"), "binary")

fplot_dat(df3a, df3a$model, or = T) 

ggsave(paste0(out_dir, "forest_mv-pooleff_bin_allmod_suppfig6A.jpg"),  width = 15, height = 8)


df3b <- sel_dat(mv.all_res, c("Pooled (FE)"), c("basic", "adjusted") , "continuous")

fplot_dat(df3b, df3b$model, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mv-pooleff_cont_allmod_suppfig6B.jpg"),  width = 15, height = 8) 




