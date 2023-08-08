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

#install.packages(c("data.table", "dplyr", "purrr", "devtools", "ggplot2"))
#devtools::install_github("NightingaleHealth/ggforestplot")

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
                mutate(
                        outname2 = ifelse(outcome == "gestational_diabetes", "Gestational diabetes", outname),
                        study = ifelse(study == "GEN3G", "Gen3G", study)
                        )

all_res <- bind_rows(mv.all_res, mr.all_res) 


############################################################################
#                      Main analyses - forest plot                         #
############################################################################

## 
## Keep only MV results for main analyses:
# a) MV adjusted models pooled across studies (fixed effects)
# b) MR IVW (not asjusted for offspring genotype) pooled across studies (fixed effects) 

main_res <- filter(all_res, study == "Pooled (FE)") %>%
            filter( 
              (analyses == "MV" & model == "adjusted") | (analyses == "MR" & model == "unadjusted" & method == "Inverse variance weighted") 
            )


############################################################################
#                           Sensitivity analyses                           #
############################################################################

##
## Directory for saving outputs

out_dir <- "FINAL_METANALYSIS_2021/output/"


##
## Function to select data subset

sel_dat <- function(dat, studies, methods, models, typevar) {

            # Filter relevant data  
            df <- filter(dat, 
                                    study %in% studies
                                    & 
                                    method %in% methods 
                                    &
                                    model %in% models 
                                    & 
                                    type %in% typevar
                           ) %>%
                    mutate(
                           method = factor(method, levels = c("MR Egger", "Weighted mode", "Weighted median", "Inverse variance weighted")),
                           group = factor(group, levels = c("Pregnancy", "Delivery", "Birth", "Postnatal"))
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


##
## Study-specific colour scheme (14 studies)

color <- c("aquamarine2", "black", "blue2", "blueviolet", "brown4", "darkcyan", 
           "chartreuse", "coral2", "darkgoldenrod3", "darkmagenta", "darkolivegreen4", "deeppink2", 
           "deepskyblue1", "burlywood3"
           )


##
## Study-specific and pooled (fixed- and random-effects) IVW results 

# Binary traits - part 1

df1a <- sel_dat(mr.all_res, unique(all_res$study), "Inverse variance weighted", "unadjusted", "binary") %>% filter(group %in% c("Pregnancy", "Delivery"))

fplot_dat(df1a, df1a$study, or = F) + xlab("logOR (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mr-studies_bin_part1.jpg"),  width = 15, height = 10)

# Binary traits - part 2

df1b <- sel_dat(mr.all_res, unique(all_res$study), "Inverse variance weighted", "unadjusted", "binary") %>% filter(group %in% c("Birth", "Postnatal"))

fplot_dat(df1b, df1b$study, or = F) + xlab("logOR (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mr-studies_bin_part2.jpg"),  width = 15, height = 10)

# Continuous traits 

df1c <- sel_dat(mr.all_res, unique(all_res$study), "Inverse variance weighted", "unadjusted", "continuous")

fplot_dat(df1c, df1c$study, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mr-studies_cont.jpg"),  width = 15, height = 10)


##
## Leave-one-out STUDY IVW results

loo_study <- function(out, vartype) {
  
  # Subset data (one outcome, all studies, only unadjusted models, one var type)
  df <- filter(mr.all_res, outcome == out
               & 
                 method == "Inverse variance weighted" 
               &
                 ! study %in% c("Pooled (RE)", "Pooled (FE)")
               & 
                 model == "unadjusted"
               &
                 type == vartype
  )
  
  if(nrow(df) == 0) {
    
    return(NULL)
    
  } else {
    
  # Format data for LOO
  study.loo_dat <- metafor::rma(yi = b, sei = se, data=df, slab = study, method="FE")
  
  # Run LOO removing one study at a time
  study.loo_res <- metafor::leave1out(study.loo_dat)
 
  # Format LOO output
  study.loo_res2 <- study.loo_res  %>% 
    as.data.frame %>%
    tibble::rownames_to_column(., var = "study") %>%
    mutate(outcome = out, se = estimate/zval )  %>%
    merge(., out_info2, by.x = "outcome", by.y = "phase1_varname") %>%
    mutate(outname2 = ifelse(outcome == "gestational_diabetes", "Gestational diabetes", outname)) %>%
    rename(b = estimate) 
  
  return(study.loo_res2)
  
  } 
}  

# Binary outcomes - part 1

loo.bin_res1 <- map(unique(mr.all_res$outcome), ~loo_study(., "binary")) %>% 
            bind_rows %>%
            filter(group %in% c("Pregnancy", "Delivery"))

fplot_dat(loo.bin_res1, loo.bin_res1$study, or = F) + xlab("logOR (95% CI)\nper 1-SD increment in BMI") + scale_color_manual(values=color)

ggsave(paste0(out_dir, "forest_loo-studies_bin_part1.jpg"),  width = 15, height = 10)

# Binary outcomes - part 2

loo.bin_res2 <- map(unique(mr.all_res$outcome), ~loo_study(., "binary")) %>% 
  bind_rows %>%
  filter(group %in% c("Birth", "Postnatal"))

fplot_dat(loo.bin_res2, loo.bin_res2$study, or = F) + xlab("logOR (95% CI)\nper 1-SD increment in BMI") + scale_color_manual(values=color)

ggsave(paste0(out_dir, "forest_loo-studies_bin_part2.jpg"),  width = 15, height = 10)

# Continuous outcomes

loo.cont_res <- map(unique(mr.all_res$outcome), ~loo_study(., "continuous")) %>% 
  bind_rows

loo.cont_res  %>%
  ggforestplot::forestplot(
    df = .,
    name = outname2,
    estimate = b,
    pvalue = pval,
    logodds = F,
    psignif = 0.05,
    xlab = "SD (95% CI)\nper 1-SD increment in BMI",
    colour = study
  ) +
  theme(
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    axis.text.x = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.y = element_text(size=15),
    strip.text = element_text(size = 17)
  )

fplot_dat(loo.cont_res, loo.cont_res$study, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI") + scale_color_manual(values=color)

ggsave(paste0(out_dir, "forest_loo-studies_cont.jpg"),  width = 15, height = 10)


##
## Multiple MR methods (pooled data only)

df2a <- sel_dat(mr.all_res, "Pooled (FE)", unique(all_res$method), "unadjusted", "binary")

fplot_dat(df2a, df2a$method, or = T) 

ggsave(paste0(out_dir, "forest_mr-methods_bin.jpg"),  width = 15, height = 8)


df2b <- sel_dat(mr.all_res, "Pooled (FE)", unique(all_res$method), "unadjusted", "continuous")

fplot_dat(df2b, df2b$method, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mr-methods_cont.jpg"),  width = 15, height = 8) 


##
## Study-specific and pooled (fixed- and random-effects) IVW results 

df3a <- sel_dat(mr.all_res, c("Pooled (FE)", "Pooled (RE)"), "Inverse variance weighted", "unadjusted", "binary")

fplot_dat(df3a, df3a$study, or = T) 

ggsave(paste0(out_dir, "forest_mr-pooleff_bin.jpg"),  width = 15, height = 8)


df3b <- sel_dat(mr.all_res, c("Pooled (FE)", "Pooled (RE)"), "Inverse variance weighted", "unadjusted", "continuous")

fplot_dat(df3b, df3b$study, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mr-pooleff_cont.jpg"),  width = 15, height = 8) 


##
## Offspring genotype unadjusted vs adjusted IVW

df4a <- filter(mr.all_res, ! group == "Pregnancy") %>%
  sel_dat(., "Pooled (FE)", "Inverse variance weighted", unique(mr.all_res$model), c("binary"))

fplot_dat(df4a, df4a$model, or = T) 

ggsave(paste0(out_dir, "forest_mr-offadj_bin.jpg"),  width = 15, height = 8)


df4b <- filter(mr.all_res, ! group == "Pregnancy") %>%
  sel_dat(., "Pooled (FE)", "Inverse variance weighted", unique(mr.all_res$model), c("continuous"))

fplot_dat(df4b, df4b$model, or = F) + xlab("SD (95% CI)\nper 1-SD increment in BMI")

ggsave(paste0(out_dir, "forest_mr-offadj_cont.jpg"),  width = 15, height = 8)


##
## SNP-specific Wald ratios

#pid <- filter(mr.all_res, study == "Pooled (FE)" & type == "binary"& model == "unadjusted") %>%
 #                 select(id.exposure, id.outcome, outcome, order) %>%
  #                distinct %>%
   #               arrange(order) %>%
    #              mutate(id = paste0(id.exposure, ".", id.outcome))%>%
     #             pull(id)

pid <- filter(mr.all_res, study == "Pooled (FE)" & model == "unadjusted") %>%
  select(id.exposure, id.outcome, outcome, type, order) %>%
  distinct %>%
  arrange(type, order) %>%
  mutate(id = paste0(id.exposure, ".", id.outcome))%>%
  pull(id)
                    
library(ggpubr)
theme_set(theme_pubr())

ggarrange(p_single[[pid[1]]], p_single[[pid[2]]], p_single[[pid[3]]], ncol = 3, nrow = 1)    %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part1.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[4]]], p_single[[pid[5]]], p_single[[pid[6]]], ncol = 3, nrow = 1)    %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part2.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[7]]], p_single[[pid[8]]], p_single[[pid[9]]], ncol = 3, nrow = 1)    %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part3.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[10]]], p_single[[pid[11]]], p_single[[pid[12]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part4.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[13]]], p_single[[pid[14]]], p_single[[pid[15]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part5.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[16]]], p_single[[pid[17]]], p_single[[pid[18]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part6.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[19]]], p_single[[pid[20]]], p_single[[pid[21]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part7.pdf"),  width = 10, height = 10)

ggarrange(p_single[[pid[22]]], p_single[[pid[23]]], p_single[[pid[24]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-wald_part8.pdf"),  width = 10, height = 10)


##
## Leave-one-out analyses

ggarrange(p_loo[[pid[1]]], p_loo[[pid[2]]], p_loo[[pid[3]]], ncol = 3, nrow = 1)    %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part1.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[4]]], p_loo[[pid[5]]], p_loo[[pid[6]]], ncol = 3, nrow = 1)    %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part2.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[7]]], p_loo[[pid[8]]], p_loo[[pid[9]]], ncol = 3, nrow = 1)    %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part3.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[10]]], p_loo[[pid[11]]], p_loo[[pid[12]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part4.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[13]]], p_loo[[pid[14]]], p_loo[[pid[15]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part5.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[16]]], p_loo[[pid[17]]], p_loo[[pid[18]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part6.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[19]]], p_loo[[pid[20]]], p_loo[[pid[21]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part7.pdf"),  width = 12, height = 10)

ggarrange(p_loo[[pid[22]]], p_loo[[pid[23]]], p_loo[[pid[24]]], ncol = 3, nrow = 1) %>% ggexport(., filename = paste0(out_dir, "forest_mr-loo_part8.pdf"),  width = 12, height = 10)


##
## Table with SNP-heterogeneity metrics (Q and MR-Egger intercept)

snp_het <- merge(out_info2, het, by.x = "phase1_varname", by.y = "outcome") %>%
            merge(., int, by = "id.outcome") %>%
            rename(egger_se = se, egger_pval = pval, exposure = exposure.x) %>%
            arrange(order) %>%
            select(exposure, outname, Q, Q_df, Q_pval, egger_intercept, egger_se, egger_pval)
  

write.csv(snp_het, file = paste0(out_dir, "snp_het_metrics.csv"))









