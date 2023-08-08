############################################################################
#                                                                          #
#                  Makes the main MV x MR & MV X PNC PLOTs                 #
#                   (and supplementary PNC plots by study)                 #
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
#                    MAIN MV X MR PLOT                                #
#---------------------------------------------------------------------#

## Load data
# MV
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")
# MR
load("./FINAL_METANALYSIS_2021/data/mr_res.RData")
##
## Combine datasets
out_info2 <- filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB))

mv.all_res <- mutate(mv.all_res, analyses = "MV") 

mr.all_res <- mutate(mr.all_res, analyses = "MR")

mr <- subset(mr.all_res , model==c("unadjusted"))
mr <- subset(mr, method==c("Inverse variance weighted"))


mr2<- mr %>% 
  group_by(outcome) %>% 
  summarise(study = study)

#remove the pooled estimates
mr2 <- subset(mr2 , study!=c("Pooled (FE)"))
mr2 <- subset(mr2 , study!=c("Pooled (RE)"))

counts <-
  mr2 %>%
  count(outcome)

names(counts)[names(counts) == 'n'] <- 'nstudies'

#merge in nstudies into the MR results
mr.all_res<-  merge(counts, mr.all_res, by = "outcome")


all_res <- bind_rows(mv.all_res, mr.all_res) %>%
  merge(., out_info2, by.x = "outcome", by.y = "phase1_varname")

#rename GDM asa character is wrong in the spelling
all_res$outname.y[all_res$outcome=="gestational_diabetes"] <- 'Gestational diabetes'

############################################################################
# Main analyses - forest plot MV X MR plot
############################################################################

# a) MV adjusted models pooled across studies (fixed effects)
# b) MR IVW (not asjusted for offspring genotype) pooled across studies (fixed effects)

main_res <- filter(all_res, study == "Pooled (FE)") %>%
  filter((analyses == "MV" & model == "adjusted") | (analyses == "MR" & model == "unadjusted" & method == "Inverse variance weighted"))


#relabel analyses for plots
main_res$analyses[main_res$analyses=="MR"] <- 'Mendelian randomization'
main_res$analyses[main_res$analyses=="MV"] <- 'Multivariable regression'

as.factor(main_res$analyses)
main_res$analyses <- factor(main_res$analyses, levels = c("Multivariable regression", 
                                                          "Mendelian randomization"))
main_res <- arrange(main_res, order.y, analyses)

main_res_b<-subset(main_res, type.y=="binary")


#Main figure - MV x MR
MVXMR_main_MA <- metagen(TE = b, seTE = se, data = main_res_b,  
                             studlab = paste(analyses), sm = "OR",
                             hakn = FALSE, byvar = outname.y, comb.fixed = F, comb.random=F)

summary(MVXMR_main_MA)

forest_MVXMR_b<-forest(MVXMR_main_MA, studlab = TRUE, 
                       comb.fixed = F,
                       type.study="square",
                       spacing = 0.7,
                       squaresize=0.5,
                       lty.fixed = 2,
                       type.fixed="square",
                       bylab = "",
                       text.fixed = "Total", # write anything 
                       text.fixed.w = "Total",
                       col.study = c("indianred", "black", "indianred", "black", "indianred", "black", 
                                     "indianred", "black", "indianred", "black", "indianred",  "black", 
                                     "indianred", "black", "indianred", "black", "indianred",  "black", 
                                     "indianred", "black", "indianred", "black", "indianred",  "black", 
                                     "indianred", "black", "indianred", "black", "indianred", "black", 
                                     "indianred", "black", "indianred", "black", "indianred", "black", 
                                     "indianred", "black", "indianred", "black" ) , 
                       col.square = c("indianred", "black", "indianred", "black", "indianred", "black", 
                                      "indianred", "black", "indianred", "black", "indianred",  "black", 
                                      "indianred", "black", "indianred", "black", "indianred",  "black", 
                                      "indianred", "black", "indianred", "black", "indianred",  "black", 
                                      "indianred", "black", "indianred", "black", "indianred", "black", 
                                      "indianred", "black", "indianred", "black", "indianred", "black", 
                                      "indianred", "black", "indianred", "black") ,
                       col.diamond="black", 
                       col.diamond.lines="black",
                       col.label.right="black",
                       col.label.left="black", 
                       colgap.right = "0.5cm",
                       colgap.left = "0.1cm",
                       colgap.forest.left ="0.1cm",
                       colgap.forest.right ="2cm",
                       col.by = "black",
                       xlab="OR per 1-SD increment in BMI", 
                       smlab=" ", 
                       leftcols=c( "studlab", "nstudies"),# To remove "logHR" and "seHR" from plot
                       leftlabs = c("Outcome", "No. of studies"),
                       rightcols=c("effect", "ci"),
                       rightlabs=c("OR","[95% CI]"),
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
                       at = c(0.5, 0.75, 1, 1.5, 2))
forest_MVXMR_b <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=450, height=5500, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/forest_MVXMR_main_b_fig2.png", forest_MVXMR_b)

############################################################################
# Main analyses - forest plot MV X PNC plot
############################################################################
#---------------------------------------------------------------------#
#                             data prep                               #
#---------------------------------------------------------------------#
##LOAD IN ALSPAC and MOBA and GENR estimates from MV (maternal main) and PNC (main)
##ALSPAC
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

##MOBA
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

##GENR
MV_main3_c <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "3 Maternal BMI CONTINUOUS")
MV_main3_c$models <-"C Maternal BMI, adjusted for paternal BMI"
MV_main3_b <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "3 Maternal BMI BINARY")
MV_main3_b$models <-"C Maternal BMI, adjusted for paternal BMI"

MV_main4_c <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "2 Paternal BMI CONTINUOUS")
MV_main4_c$models <-"D Paternal BMI, adjusted for maternal BMI"
MV_main4_b <- read_excel("./GENERATIONR/data_analysis/results_files/Maternal and paternal BMI analysis bio.xls" , sheet = "2 Paternal BMI BINARY")
MV_main4_b$models <-"D Paternal BMI, adjusted for maternal BMI"

#put together

MV_main_b_GENR <- bind_rows(MV_main3_b, MV_main4_b)
MV_main_c_GENR <- bind_rows(MV_main3_c, MV_main4_c)
MV_main_b_GENR$abbrev <-"GENR"
MV_main_c_GENR$abbrev <-"GENR"

MV_main_b <- bind_rows(MV_main_b_MOBA, MV_main_b_ALSPAC, MV_main_b_GENR)
MV_main_c <- bind_rows(MV_main_c_MOBA, MV_main_c_ALSPAC, MV_main_c_GENR)

#---------------------------------------------------------------------#
#               merge PNC & MV data with outcomes names               #
#---------------------------------------------------------------------#
load("./FINAL_METANALYSIS_2021/data/mv_res.RData")
out_info2 <- filter(out_info, primary == "yes") %>%
  select(!(phase1_units:phase2_varname_UKBB))

MV_main_b1<-MV_main_b
names(MV_main_b1)[names(MV_main_b1) == 'outcome'] <- 'phase2_varname'
res_MV_b1 <- merge(MV_main_b1, out_info, by = "phase2_varname")
names(res_MV_b1)[names(res_MV_b1) == 'phase2_varname'] <- 'varname'
res_MV_b1 <- subset(res_MV_b1 , select=c("group", "outname", "models", "abbrev", "N", "cases", "logOR", "se", "lci", "uci", "pval", "type", "varname", "order"))
res_MV_b1$outname[res_MV_b1$varname=="gdm_subsamp"] <- 'Gestational diabetes'

MV_main_b2<-MV_main_b
names(MV_main_b2)[names(MV_main_b2) == 'outcome'] <- 'phase1_varname'
res_MV_b2 <- merge(MV_main_b2, out_info2, by = "phase1_varname") 
names(res_MV_b2)[names(res_MV_b2) == 'phase1_varname'] <- 'varname'
res_MV_b2 <- subset(res_MV_b2 , select=c("group", "outname", "models", "abbrev", "N", "cases", "logOR", "se", "lci", "uci", "pval", "type", "varname", "order"))
res_MV_b2$outname[res_MV_b2$varname=="gdm_subsamp"] <- 'Gestational diabetes'


res_MV_b<-rbind(res_MV_b1, res_MV_b2)
#---------------------------------------------------------------------#
#                           PNC by study                              #
#---------------------------------------------------------------------#

res_MV<-res_MV_b
names(res_MV)[names(res_MV) == 'logOR'] <- 'b'

as.factor(res_MV$group)

res_MV <- arrange(res_MV, order, abbrev, group)

res_MV_mat <- subset(res_MV, models=="C Maternal BMI, adjusted for paternal BMI")
res_MV_pnc <- subset(res_MV , models!="C Maternal BMI, adjusted for paternal BMI")


MV_main_MA <- metagen(TE = b, seTE = se, data = res_MV_mat, 
                      studlab = paste(abbrev), sm = "OR",
                      hakn = FALSE, byvar = outname, comb.fixed = T, comb.random=F)


print(MV_main_MA)

#forest plot - maternal BMI adjusted for pat BMI
forest_mat_b<-forest(MV_main_MA, studlab = TRUE, 
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
                     leftcols=c("studlab", "cases", "N"),# To remove "logHR" and "seHR" from plot
                     leftlabs = c("Study", "Cases", "N"),
                     rightcols=c("effect", "ci", "w.fixed"),
                     rightlabs=c("OR","[95% CI]", "Weight"),
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
                     at = c(0.5, 0.75, 1, 1.5, 2))
forest_mat_b <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=6000, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/studyspecificPNC_mat_b_supp12A.png", forest_mat_b)



#forest plot - paternal BMI adjusted for mat BMI by study
PNC_main_MA <- metagen(TE = b, seTE = se, data = res_MV_pnc, 
                       studlab = paste(abbrev), sm = "OR",
                       hakn = FALSE, byvar = outname, comb.fixed = T, comb.random=F)


print(PNC_main_MA)


#forest plot - maternal
forest_pnc_b<-forest(PNC_main_MA, studlab = TRUE, 
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
                     leftcols=c("studlab", "cases", "N"),# To remove "logHR" and "seHR" from plot
                     leftlabs = c("Study", "Cases", "N"),
                     rightcols=c("effect", "ci", "w.fixed"),
                     rightlabs=c("OR","[95% CI]", "Weight"),
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
                     at = c(0.5, 0.75, 1, 1.5, 2))
forest_pnc_b <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=6000, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/studyspecificPNC_pnc_b_supp13A.png", forest_pnc_b)

#---------------------------------------------------------------------#
#                    MAIN MV X PNC PLOT                                #
#---------------------------------------------------------------------#

#add total cases and total N
res_MV_mat_total<- res_MV_mat %>% 
  group_by(outname) %>% 
  summarise(Totalcases = sum(cases), TotalN = sum(N), outname=outname,  order=order, group=group, abbrev=abbrev)

res_MV_mat_total<- res_MV_mat_total %>% 
  add_count(outname) 

#, order=order, group=group, abbrev=abbrev
res_MV_pnc_total<- res_MV_pnc %>% 
  group_by(outname) %>% 
  summarise(Totalcases = sum(cases), TotalN = sum(N), outname=outname,  order=order, group=group, abbrev=abbrev)


res_MV_pnc_total<- res_MV_pnc_total %>% 
  add_count(outname) 


#pooled data
pooled.mv <- data.frame(summary(MV_main_MA)$bylevs, summary(MV_main_MA)$TE.fixed.w, summary(MV_main_MA)$seTE.fixed.w)
names(pooled.mv)[1] <- "Outcome"
names(pooled.mv)[2] <- "lnor"
names(pooled.mv)[3] <- "se"
pooled.mv$group <- "Multivariable regression"

pooled.pnc <- data.frame(summary(PNC_main_MA)$bylevs, summary(PNC_main_MA)$TE.fixed.w, summary(PNC_main_MA)$seTE.fixed.w)
names(pooled.pnc)[1] <- "Outcome"
names(pooled.pnc)[2] <- "lnor"
names(pooled.pnc)[3] <- "se"
pooled.pnc$group <- "Paternal negative control"

res_MV_mat_total<-subset(res_MV_mat_total, abbrev=="ALSPAC")
names(res_MV_mat_total)[names(res_MV_mat_total) == 'outname'] <- 'Outcome'
res_MV_mat_total<-subset(res_MV_mat_total, select=c("Outcome", "Totalcases", "TotalN", "order", "group", "n"))
res.mv.mat <- merge(pooled.mv, res_MV_mat_total, by = "Outcome")

res_MV_pnc_total<-subset(res_MV_pnc_total, abbrev=="ALSPAC")
names(res_MV_pnc_total)[names(res_MV_pnc_total) == 'outname'] <- 'Outcome'
res_MV_pnc_total<-subset(res_MV_pnc_total, select=c("Outcome", "Totalcases", "TotalN", "order", "group", "n"))
res.mv.pnc <- merge(pooled.pnc, res_MV_pnc_total, by = "Outcome")

pooled.mv.pnc<-rbind(res.mv.mat,res.mv.pnc)

as.factor(pooled.mv.pnc$group.x)

pooled.mv.pnc <- arrange(pooled.mv.pnc, order, group.x)

pooled.mv.pnc$Outcome[pooled.mv.pnc$phase2_varname=="Gestational? diabetes"] <- 'Gestational diabetes'


MVXPNC_main_MA_np <- metagen(TE = lnor, seTE = se, data = pooled.mv.pnc, 
                             studlab = paste(group.x), sm = "OR",
                             hakn = FALSE, byvar = Outcome, comb.fixed = F, comb.random=F)
print(MVXPNC_main_MA_np)

#forest plot
forest_MVXPNC_b_np<-forest(MVXPNC_main_MA_np, studlab = TRUE, 
                           comb.fixed = F,
                           type.study="square",
                           spacing = 0.80,
                           squaresize=0.5,
                           lty.fixed = 2,
                           type.fixed="square",
                           bylab = "",
                           text.fixed = "Total", # write anything 
                           text.fixed.w = "Total",
                           col.study = c("indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                         "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue", 
                                         "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                         "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue",
                                         "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                         "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue",
                                         "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                         "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue") , 
                           col.square = c("indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                          "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue", 
                                          "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                          "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue",
                                          "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                          "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue",
                                          "indianred", "deepskyblue", "indianred", "deepskyblue", "indianred",
                                          "deepskyblue", "indianred", "deepskyblue", "indianred", "deepskyblue") ,
                           col.diamond="black", 
                           col.diamond.lines="black",
                           col.label.right="black",
                           col.label.left="black", 
                           colgap.right = "0.5cm",
                           colgap.left = "0.1cm",
                           colgap.forest.left ="0.1cm",
                           colgap.forest.right ="2cm",
                           col.by = "black",
                           xlab="OR per 1-SD increment in BMI", 
                           smlab=" ", 
                           leftcols=c( "studlab", "n"),# To remove "logHR" and "seHR" from plot
                           leftlabs = c("Outcome", "No. of studies"),
                           rightcols=c("effect", "ci"),
                           rightlabs=c("OR","[95% CI]"),
                           just = "c", just.addcols = c("c"), just.studlab = "l",
                           test.overall = F,
                           lwd=2,
                           print.I2 = F,
                           plotwidth="8.0cm",
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
                           at = c(0.50, 0.75, 1, 1.25, 1.5, 2))
forest_MVXPNC_b_np <- recordPlot()

#save
save_func <- function(file_name, plot_name)
{
  png(file_name, res=400, height=6000, width=3500)
  print(plot_name)
  dev.off()
  
}

save_func("FINAL_METANALYSIS_2021/output/forest_MVXPNC_main_b_fig3.png", forest_MVXPNC_b_np)

