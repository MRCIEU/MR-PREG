############################################################################
#                                                                          #
#                 This script creates supplemtary tables                   #
#                                                                          #
############################################################################
#---------------------------------------------------------------------#
#                            Housekeeping                             #
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 
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

#MAIN MV X MR PLOT
setwd("./")

############################################################################
# Format data #
############################################################################

## Load data
# MV
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")
# MR
load("./FINAL_METANALYSIS_2021/data/mr_res.RData")

# MV
mv <- subset(mv.all_res , study==c("Pooled (FE)", "Pooled (RE)"))

mv <- subset(mv , select=c("outcome", "outname", "study", "model", "order", "group", "q", "q_pval", "i2", "nstudies", "type"))

mv <- subset(mv, study==c("Pooled (FE)"))
mv <- subset(mv, model==c("adjusted"))


mv <- arrange(mv, type, order)

mv_het <- subset(mv , select=c("group", "outname", "nstudies", "q", "q_pval", "i2"))

write.xlsx(mv_het, file = "FINAL_METANALYSIS_2021/output/mv_het.xlsx", overwrite=T)
mv_het <- read.xlsx(xlsxFile = "FINAL_METANALYSIS_2021/output/mv_het.xlsx")



