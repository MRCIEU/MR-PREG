############################################################################
#                                                                          #
#          This script pools MV association estimates across studies       #
#                    split by pre/during pregnancy weight                  #
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
library(haven)
library(meta)


############################################################################
#                          Import auxiliary files                          #
############################################################################

##
## Study info file

studies <- read.csv("./FINAL_METANALYSIS_2021/data/list_of_studies.csv")

studies.p1 <- filter(studies, phase1 == 1 & keep == "yes") %>% pull(study_var)


##
## Outcomes info file 

out_info <- read.csv("./FINAL_METANALYSIS_2021/data/list_of_varnames.csv") %>%
  rename(outname = "?..Outcome") %>%
  filter(primary == "yes")


############################################################################
#                   Import and format files - phase 1                      #
############################################################################

##
## Import MV estimates for phase 1 studies

mv.p1_extr <- function(study, covar) {
  
  df <- read_dta(paste0("./META_ANALYSIS/single_study_results/", study, "_obs_", covar, ".dta")) 
  
  if("p" %in% names(df)) {
    
    df <- rename(df, pval = p)
    
  }
  
  df <- df %>%
    rename(study = Study, b = beta) %>%
    mutate(model = covar) %>%
    merge(., out_info, by.x = "outcome", by.y = "phase1_varname") %>%
    select(outcome, b, se, pval, study, model, N, outname, primary, type, order, group) 
}           

mv.p1_args <-  cross_df(list( 
  study = studies.p1[ ! studies.p1 %in% c("ALSPAC", "UK_BIOBANK") ] , 
  covar = c("basic", "adjusted") 
))

mv.p1_res <- pmap(mv.p1_args, ~mv.p1_extr(.x, .y)) %>% bind_rows

# P.S.: Two studies were removed: ALSPAC (also present in phase 2) and UKB (no data on maternal BMI)

############################################################################
#                   Import and format files - phase 2                      #
############################################################################

##
## Import MV estimates for phase 2 studies

mv.p2_extr <- function(study, vartype, adj) {
  
  # Import MV data 
  if(study == "ALSPAC") {
    
    fpath <- paste0(study, "/results_files/Maternal and paternal BMI analysis.xls")
    
  } else {
    
    fpath <- paste0(study, "/results/Maternal and paternal BMI analysis.xls")
  }
  
  # Merge with info file
  df <- readxl::read_excel(fpath, sheet = paste0(vartype, " ", adj)) %>%
    mutate(study = study, model = tolower(adj)) %>%
    merge(., out_info, by.x = "outcome", by.y = "phase2_varname") %>%
    mutate(outcome = phase1_varname)
  
  # Harmonise variable names
  if("beta" %in% names(df)) { 
    
    df <- rename(df, b = beta)
    
  } else if("logOR" %in% names(df)) {
    
    df <- rename(df, b = logOR)
    
  }
  
  # Select variables
  df <- df %>%
    select(outcome, b, se, pval, study, model, N, outname, primary, type, order, group) 
  
}

mv.p2_args <-  cross_df(list( 
  study = c("ALSPAC", "MOBA") , 
  vartype = c("Binary", "Continuous"),
  adj = c("Basic", "Adjusted")
))

mv.p2_res <- pmap(mv.p2_args, mv.p2_extr) %>% bind_rows

# P.S.: Only ALSPAC and MoBa had maternal BMI data among phase 2 studies

############################################################################
#                   Metanalyse MV estimates across studies                 #
############################################################################

##
## Combine MV data from phase 1 and 2 studies

mv_res <- bind_rows(mv.p1_res, mv.p2_res)

#Generate an indicator variable
mv_res <- mv_res %>% mutate(mv_res, ind = case_when(study == "BIB" | study == "GENR" |study == "HAPO" | study == "MOBA" ~ 1,
      study == "GOYA" | study == "ALSPAC" |study == "DNBC" | study == "EFSOCH"| study == "GEN3G" |study == "INMA" | study == "NFBC1966" | study == "NFBC1986" ~ 0))
      
#change NAs to 0
mv_res$ind[is.na(mv_res$ind)] <- 0


##
## Function to run metanalyses

mv_res_1<-subset(mv_res, ind==1)
mv_res_1$ind<-"During pregnancy"

mv_res_0<-subset(mv_res, ind==0)
mv_res_0$ind<-"Pre-pregnancy"

mv_res <- bind_rows(mv_res_1, mv_res_0)


#ind==1
mv_meta_1 <- function(out, adj) {
  
  # Keep one outcome and one model at a time
  input <- filter(mv_res_1, outcome == out & model == adj) 
  print(input)
  
  # Meta-analysis across studies
  ma <- metagen(TE = b, 
                seTE = se, 
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
    b             = ma$TE.fixed,
    se            = ma$seTE.fixed,
    pval          = ma$pval.fixed,
    nstudies      = ma$k,
    q             = ma$Q,
    q_pval        = ma$pval.Q,
    i2            = ma$I2,
    study         = "Pooled (FE)"
  )
  
  # Extract values from random-effect metanalysis
  ran_ma <- data.table(
    b            = ma$TE.random,
    se           = ma$seTE.random,
    pval         = ma$pval.random,
    nstudies     = ma$k,
    q            = ma$Q,
    q_pval       = ma$pval.Q,
    i2           = ma$I2,
    study        = "Pooled (RE)"
  )
  
  # Combine results fom fixed- and random-effects metanalyses
  ma <- bind_rows(fix_ma, ran_ma) %>%
        mutate(exposure = "Maternal BMI", outcome = out, model = adj)
  
  return(ma)
}


mv.pool_args <-  cross_df(list( 
  out = unique(mv_res_1$outcome) , 
  adj = c("basic", "adjusted")
))



mv.pool_res_1 <- pmap(mv.pool_args, mv_meta_1) %>% 
  bind_rows %>%
  merge(., out_info, by.x = "outcome", by.y = "phase1_varname") %>%
  select(outcome, b, se, pval, q, q_pval, i2, study, model, nstudies, outname, primary, type, order, group) 
   


#ind==0
mv_meta_0 <- function(out, adj) {
  
  # Keep one outcome and one model at a time
  input <- filter(mv_res_0, outcome == out & model == adj) 
  print(input)
  
  # Meta-analysis across studies
  ma <- metagen(TE = b, 
                seTE = se, 
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
    b             = ma$TE.fixed,
    se            = ma$seTE.fixed,
    pval          = ma$pval.fixed,
    nstudies      = ma$k,
    q             = ma$Q,
    q_pval        = ma$pval.Q,
    i2            = ma$I2,
    study         = "Pooled (FE)"
  )
  
  # Extract values from random-effect metanalysis
  ran_ma <- data.table(
    b            = ma$TE.random,
    se           = ma$seTE.random,
    pval         = ma$pval.random,
    nstudies     = ma$k,
    q            = ma$Q,
    q_pval       = ma$pval.Q,
    i2           = ma$I2,
    study        = "Pooled (RE)"
  )
  
  # Combine results fom fixed- and random-effects metanalyses
  ma <- bind_rows(fix_ma, ran_ma) %>%
    mutate(exposure = "Maternal BMI", outcome = out, model = adj)
  
  return(ma)
}


mv.pool_args <-  cross_df(list( 
  out = unique(mv_res_0$outcome) , 
  adj = c("basic", "adjusted")
))



mv.pool_res_0 <- pmap(mv.pool_args, mv_meta_0) %>% 
  bind_rows %>%
  merge(., out_info, by.x = "outcome", by.y = "phase1_varname") %>%
  select(outcome, b, se, pval, q, q_pval, i2, study, model, nstudies, outname, primary, type, order, group) 

##combine the two preg weight results
mv.pool_res_1$ind<-"During pregnancy"
mv.pool_res_0$ind<-"Pre-pregnancy"
mv.pool_res <- bind_rows(mv.pool_res_1, mv.pool_res_0)

#c("Pre-pregnancy", "During pregnancy")

## Combine all results

mv.all_res <- bind_rows(mv_res, mv.pool_res)


##
## Save relevant objects

save(list = c("studies", "out_info", "mv.all_res"), file = "./FINAL_METANALYSIS_2021/data/mv_res_wgt.RData")

write.csv(mv.all_res, file = "./FINAL_METANALYSIS_2021/data/mv_all_res_wgt.csv")




