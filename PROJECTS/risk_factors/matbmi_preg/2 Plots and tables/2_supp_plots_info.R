############################################################################
#                                                                          #
#        Makes the supplementary main MV x MR & MV X PNC PLOTs             #
#                         for continuous outcomes                          #
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

setwd("./")

#---------------------------------------------------------------------#
#                        MAIN MV X MR PLOT                            #
#                       continuous outcomes                           #
#---------------------------------------------------------------------#

## Load data
# MV
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")
# MR
load("./FINAL_METANALYSIS_2021/data/mr_res.RData")

## Combine datasets
out_info2 <- filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB))

mv.all_res <- mutate(mv.all_res, analyses = "MV") 

mr.all_res <- mutate(mr.all_res, analyses = "MR")

#add in Ns for mR studies
mr <- subset(mr.all_res , model==c("unadjusted"))
mr <- subset(mr, method==c("Inverse variance weighted"))


mr2<- mr %>% 
  group_by(outcome) %>% 
  summarise(study = study)


mr2 <- subset(mr2 , study!=c("Pooled (FE)"))
mr2 <- subset(mr2 , study!=c("Pooled (RE)"))

counts <-
  mr2 %>%
  count(outcome)

names(counts)[names(counts) == 'n'] <- 'nstudies'

#
mr.all_res<-  merge(counts, mr.all_res, by = "outcome")



all_res <- bind_rows(mv.all_res, mr.all_res) %>%
  merge(., out_info2, by.x = "outcome", by.y = "phase1_varname")

############################################################################
# Main analyses - forest plot #
############################################################################

# a) MV adjusted models pooled across studies (fixed effects)
# b) MR IVW (not asjusted for offspring genotype) pooled across studies (fixed effects)

main_res <- filter(all_res, study == "Pooled (FE)") %>%
  filter((analyses == "MV" & model == "adjusted") | (analyses == "MR" & model == "unadjusted" & method == "Inverse variance weighted"))

main_res$analyses[main_res$analyses=="MR"] <- 'Mendelian randomization'
main_res$analyses[main_res$analyses=="MV"] <- 'Multivariable regression'

as.factor(main_res$analyses)
main_res$analyses <- factor(main_res$analyses, levels = c("Multivariable regression", 
                                                          "Mendelian randomization"))
main_res <- arrange(main_res, order.y, analyses)

main_res_c<-subset(main_res, type.y=="continuous")


MVXMR_main_MA <- metagen(TE = b, seTE = se, data = main_res_c, 
                         studlab = paste(outname.y), sm = "MD",
                         hakn = FALSE, byvar = group.y, comb.fixed = F, comb.random=F)

MVXMR_main_MA <- metagen(TE = b, seTE = se, data = main_res_c,  
                             studlab = paste(analyses), sm = "MD",
                             hakn = FALSE, byvar = outname.y, comb.fixed = F, comb.random=F)

print(MVXMR_main_MA)

#forest plot - maternal
forest_MVXMR_c<-forest(MVXMR_main_MA, studlab = TRUE, 
                       comb.fixed = F,
                       type.study="square",
                       spacing = 0.7,
                       squaresize=0.5,
                       lty.fixed = 2,
                       type.fixed="square",
                       bylab = "",
                       text.fixed = "Total", # write anything 
                       text.fixed.w = "Total",
                       col.study = c("indianred", "black", "indianred", "black", 
                                     "indianred",  "black", "indianred", "black" ) , 
                       col.square = c("indianred", "black", "indianred", "black", 
                                      "indianred", "black", "indianred", "black" ) ,
                       col.diamond="black", 
                       col.diamond.lines="black",
                       col.label.right="black",
                       col.label.left="black", 
                       colgap.right = "0.5cm",
                       colgap.left = "0.1cm",
                       colgap.forest.left ="0.1cm",
                       colgap.forest.right ="2cm",
                       col.by = "black",
                       xlab="MD per 1-SD increment in BMI", 
                       smlab=" ", 
                       leftcols=c( "studlab", "nstudies"),# To remove "logHR" and "seHR" from plot
                       leftlabs = c("Outcome", "No. of studies"),
                       rightcols=c("effect", "ci"),
                       rightlabs=c("MD","[95% CI]"),
                       just = "c", just.addcols = c("c"), just.studlab = "l",
                       test.overall = F,
                       lwd=2,
                       print.I2 = F,
                       plotwidth="5.5cm",
                       print.I2.ci = FALSE, 
                       print.tau2 = F, 
                       print.Q = FALSE,
                       digits = 2, 
                       fontsize = 11,
                       overall = FALSE,
                       overall.hetstat = FALSE,
                       print.subgroup.name=F,
                       hetstat=F,
                       test.subgroup.fixed=FALSE,
                       at = c(-0.1, 0, 0.1, 0.2, 0.3))
forest_MVXMR_c <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=4500, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/forest_MVXMR_main_c_suppfig1.png", forest_MVXMR_c)
#---------------------------------------------------------------------#
#                             PNC data                                #
#---------------------------------------------------------------------#

##LOAD IN ALSPAC
MV_main3_c <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "3 Maternal BMI CONTINUOUS")
MV_main3_c$models <-"C Maternal BMI, adjusted for paternal BMI"
MV_main3_b <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "3 Maternal BMI BINARY")
MV_main3_b$models <-"C Maternal BMI, adjusted for paternal BMI"

MV_main4_c <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "2 Paternal BMI CONTINUOUS")
MV_main4_c$models <-"D Paternal BMI, adjusted for maternal BMI"
MV_main4_b <- read_excel("./ALSPAC/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "2 Paternal BMI BINARY")
MV_main4_b$models <-"D Paternal BMI, adjusted for maternal BMI"

MV_main_b_ALSPAC <- bind_rows(MV_main3_b, MV_main4_b)
MV_main_c_ALSPAC <- bind_rows(MV_main3_c, MV_main4_c)
MV_main_b_ALSPAC$abbrev <-"ALSPAC"
MV_main_c_ALSPAC$abbrev <-"ALSPAC"

##LOAD IN MOBA
MV_main3_c <- read_excel("./MOBA/results/Maternal and paternal BMI analysis.xls" , sheet = "3 Maternal BMI CONTINUOUS")
MV_main3_c$models <-"C Maternal BMI, adjusted for paternal BMI"
MV_main3_b <- read_excel("./MOBA/results/Maternal and paternal BMI analysis.xls" , sheet = "3 Maternal BMI BINARY")
MV_main3_b$models <-"C Maternal BMI, adjusted for paternal BMI"

MV_main4_c <- read_excel("./MOBA/results/Maternal and paternal BMI analysis.xls" , sheet = "2 Paternal BMI CONTINUOUS")
MV_main4_c$models <-"D Paternal BMI, adjusted for maternal BMI"
MV_main4_b <- read_excel("./MOBA/results/Maternal and paternal BMI analysis.xls" , sheet = "2 Paternal BMI BINARY")
MV_main4_b$models <-"D Paternal BMI, adjusted for maternal BMI"

#put mv_main together

MV_main_b_MOBA <- bind_rows(MV_main3_b, MV_main4_b)
MV_main_c_MOBA <- bind_rows(MV_main3_c, MV_main4_c)
MV_main_b_MOBA$abbrev <-"MOBA"
MV_main_c_MOBA$abbrev <-"MOBA"

##LOAD IN GENR
MV_main3_c <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "3 Maternal BMI CONTINUOUS")
MV_main3_c$models <-"C Maternal BMI, adjusted for paternal BMI"
MV_main3_b <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "3 Maternal BMI BINARY")
MV_main3_b$models <-"C Maternal BMI, adjusted for paternal BMI"

MV_main4_c <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "2 Paternal BMI CONTINUOUS")
MV_main4_c$models <-"D Paternal BMI, adjusted for maternal BMI"
MV_main4_b <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "2 Paternal BMI BINARY")
MV_main4_b$models <-"D Paternal BMI, adjusted for maternal BMI"

#put mv_main together

MV_main_b_GENR <- bind_rows(MV_main3_b, MV_main4_b)
MV_main_c_GENR <- bind_rows(MV_main3_c, MV_main4_c)
MV_main_b_GENR$abbrev <-"GENR"
MV_main_c_GENR$abbrev <-"GENR"

MV_main_b <- bind_rows(MV_main_b_MOBA, MV_main_b_ALSPAC, MV_main_b_GENR)
MV_main_c <- bind_rows(MV_main_c_MOBA, MV_main_c_ALSPAC, MV_main_c_GENR)


#merge with Carolinas dataset for outcome vars
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")
out_info2 <- filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB))

mv.all_res <- mutate(mv.all_res, analyses = "MV") 
#BINARY
#merge with Carolinas outcomes

names(MV_main_b)[names(MV_main_b) == 'outcome'] <- 'phase2_varname'
res_MV_b <- merge(MV_main_b, out_info, by = "phase2_varname")

#CONTINUOUS
#merge with Carolinas outcomes

MV_main_c$outcome[MV_main_c$outcome=="ponderal_zscore"] <- 'ponderal_index'
MV_main_c$outcome[MV_main_c$outcome=="birthlength_zscore"] <- 'birthlength'
MV_main_c$outcome[MV_main_c$outcome=="zbw_subsamp"] <- 'birthweight'
MV_main_c$outcome[MV_main_c$outcome=="birthweight_zscore"] <- 'birthweight'
MV_main_c$outcome[MV_main_c$outcome=="gestage_zscore"] <- 'gestational_age'

names(MV_main_c)[names(MV_main_c) == 'outcome'] <- 'phase1_varname'

res_MV_c <- merge(MV_main_c, out_info, by = "phase1_varname") 

as.factor(res_MV_c$abbrev)

res_MV_c <- arrange(res_MV_c, order, abbrev, group)


res_MV_main_c <- subset(res_MV_c , models=="C Maternal BMI, adjusted for paternal BMI")
res_MV_pnc_c <- subset(res_MV_c , models!="C Maternal BMI, adjusted for paternal BMI")

#---------------------------------------------------------------------#
#               meta-analyse mat and PNC                              #
#---------------------------------------------------------------------#
res_MV_main_c <- arrange(res_MV_main_c, order, abbrev)

MV_main_MA <- metagen(TE = beta, seTE = se, data = res_MV_main_c, 
                      studlab = paste(abbrev), sm = "MD",
                      hakn = FALSE, byvar = outname, comb.fixed = T, comb.random=F)


print(MV_main_MA)


#forest plot - maternal
forest_mat_c<-forest(MV_main_MA, studlab = TRUE, 
                     comb.fixed = MV_main_MA$comb.fixed,
                     type.study="square",
                     spacing = 0.62,
                     squaresize=0.5,
                     lty.fixed = 2,
                     type.fixed="square",
                     bylab = "",
                     text.fixed = "Total", # write anything 
                     text.fixed.w = "Total",
                     col.study="black", 
                     col.square="black", 
                     col.diamond="indianred", 
                     col.diamond.lines="indianred",
                     col.label.right="black",
                     col.label.left="black", 
                     colgap.right = "0.5cm",
                     colgap.forest.left ="2.2cm",
                     col.by = "black",
                     smlab="", 
                     leftcols=c("studlab", "N"),# To remove "logHR" and "seHR" from plot
                     leftlabs = c("Study", "N"),
                     rightcols=c("effect", "ci", "w.fixed"),
                     rightlabs=c("MD","[95% CI]", "Weight"),
                     test.overall = F,
                     lwd=0.1,
                     print.I2 = F,
                     plotwidth="7.5cm",
                     print.I2.ci = FALSE, 
                     print.tau2 = F, 
                     print.Q = FALSE,
                     digits = 2, 
                     fontsize = 9,
                     overall = FALSE,
                     overall.hetstat = FALSE,
                     print.subgroup.name=F,
                     hetstat=F,
                     test.subgroup.fixed=FALSE,
                     at = c(-0.2, -0.1, 0, 0.1, 0.2))
forest_mat_c <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=6000, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/studyspecificPNC_mat_c_supp12B.png", forest_mat_c)



#PNC
res_MV_pnc_c <- arrange(res_MV_pnc_c, order, abbrev)

PNC_main_MA <- metagen(TE = beta, seTE = se, data = res_MV_pnc_c, 
                       studlab = paste(abbrev), sm = "MD",
                       hakn = FALSE, byvar = outname, comb.fixed = T, comb.random=F)


print(PNC_main_MA)


#forest plot - maternal
forest_pnc_c<-forest(PNC_main_MA, studlab = TRUE, 
                     comb.fixed = PNC_main_MA$comb.fixed,
                     type.study="square",
                     spacing = 0.62,
                     squaresize=0.5,
                     lty.fixed = 2,
                     type.fixed="square",
                     bylab = "",
                     text.fixed = "Total", # write anything 
                     text.fixed.w = "Total",
                     col.study="black", 
                     col.square="black", 
                     col.diamond="deepskyblue", 
                     col.diamond.lines="deepskyblue",
                     col.label.right="black",
                     col.label.left="black", 
                     colgap.right = "0.5cm",
                     colgap.forest.left ="2.2cm",
                     col.by = "black",
                     smlab="", 
                     leftcols=c("studlab", "N"),# To remove "logHR" and "seHR" from plot
                     leftlabs = c("Study", "N"),
                     rightcols=c("effect", "ci", "w.fixed"),
                     rightlabs=c("MD","[95% CI]", "Weight"),
                     test.overall = F,
                     lwd=0.1,
                     print.I2 = F,
                     plotwidth="7.5cm",
                     print.I2.ci = FALSE, 
                     print.tau2 = F, 
                     print.Q = FALSE,
                     digits = 2, 
                     fontsize = 9,
                     overall = FALSE,
                     overall.hetstat = FALSE,
                     print.subgroup.name=F,
                     hetstat=F,
                     test.subgroup.fixed=FALSE,
                     at = c( -0.05, 0, 0.05))
forest_pnc_c <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=6000, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/studyspecificPNC_pnc_c_supp13B.png", forest_pnc_c)

#---------------------------------------------------------------------#
#                    MAIN MV X PNC PLOT                               #
#---------------------------------------------------------------------#

#pooled data
pooled.mv <- data.frame(summary(MV_main_MA)$bylevs, summary(MV_main_MA)$TE.fixed.w, summary(MV_main_MA)$seTE.fixed.w)
names(pooled.mv)[1] <- "Outcome"
names(pooled.mv)[2] <- "b"
names(pooled.mv)[3] <- "se"
pooled.mv$group <- "Multivariable regression"

pooled.pnc <- data.frame(summary(PNC_main_MA)$bylevs, summary(PNC_main_MA)$TE.fixed.w, summary(PNC_main_MA)$seTE.fixed.w)
names(pooled.pnc)[1] <- "Outcome"
names(pooled.pnc)[2] <- "b"
names(pooled.pnc)[3] <- "se"
pooled.pnc$group <- "Paternal negative control"

pooled.mv.pnc<-rbind(pooled.mv,pooled.pnc)

#merge with Carolinas dataset
main_res2 <- filter(all_res, study == "Pooled (FE)") %>%
  filter((analyses == "MV" & model == "adjusted"))
  
names(main_res2)[names(main_res2) == 'outname.y'] <- 'Outcome'


res <- merge(pooled.mv.pnc, main_res2, by = "Outcome") 


as.factor(res$group)

res <- arrange(res, order.y, group)

res$nstudies2<-"3"

MVXPNC_main_MA_np <- metagen(TE = b.x, seTE = se.x, data = res, 
                             studlab = paste(group), sm = "MD",
                             hakn = FALSE, byvar = Outcome, comb.fixed = F, comb.random=F)


print(MVXPNC_main_MA_np)

#forest plot - maternal
forest_MVXPNC_c_np<-forest(MVXPNC_main_MA_np, studlab = TRUE, 
                           comb.fixed = F,
                           type.study="square",
                           spacing = 0.9,
                           squaresize=0.5,
                           lty.fixed = 2,
                           type.fixed="square",
                           bylab = "",
                           text.fixed = "Total", # write anything 
                           text.fixed.w = "Total",
                           col.study = c("indianred", "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue") , 
                           col.square = c("indianred", "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue") ,
                           col.diamond="black", 
                           col.diamond.lines="black",
                           col.label.right="black",
                           col.label.left="black", 
                           colgap.right = "0.5cm",
                           colgap.left = "0.1cm",
                           colgap.forest.left ="0.1cm",
                           colgap.forest.right ="2cm",
                           col.by = "black",
                           xlab="MD per 1-SD increment in BMI", 
                           smlab=" ", 
                           leftcols=c( "studlab", "nstudies2"),# To remove "logHR" and "seHR" from plot
                           leftlabs = c("Outcome", "No. of studies"),
                           rightcols=c("effect", "ci"),
                           rightlabs=c("MD","[95% CI]"),
                           just = "c", just.addcols = c("c"), just.studlab = "l",
                           test.overall = F,
                           lwd=2,
                           print.I2 = F,
                           plotwidth="7.5cm",
                           print.I2.ci = FALSE, 
                           print.tau2 = F, 
                           print.Q = FALSE,
                           digits = 2, 
                           fontsize = 12,
                           overall = FALSE,
                           overall.hetstat = FALSE,
                           print.subgroup.name=F,
                           hetstat=F,
                           test.subgroup.fixed=FALSE,
                           at = c(-0.1, -0.05, 0, 0.05, 0.1))
forest_MVXPNC_c_np <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=4500, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/forest_MVXPNC_main_c_suppfig3.png", forest_MVXPNC_c_np)

