# Clear the work environment
rm(list = ls())

# Load packages
require(tidyverse)
require(data.table)

####################################################################################
##                                       Input                                    ##
####################################################################################

# Directories
sumdat_dir=paste0(Sys.getenv('MRPREG_sumdat'), 
  "results/GWAS/")

# List of MR-PREG outcomes 
varlist_fpath <- paste0(sumdat_dir, "Metanalysis/mothers/", "vars_metanalysis_v5_mum_final.txt")
varlist <- read.table(varlist_fpath, header = T, sep = "\t")
# P.S.: include only studies eligible for GWAMA

# Split into binary and continuous outcomes
cont.out <- c("zbw_all","apgar1","apgar5","bf_dur_4c","nvp_sev_all","nvp_sev_subsamp","ga_all", "ga_subsamp")
bin.out <- filter(varlist, !Variable_name %in% cont.out) %>% select(Variable_name) %>% pull

# Variables to Extract
cont.vars <- c("SNP", "samplesize")
bin.vars <- c("SNP", "ncase", "ncontrol", "samplesize")
meta.vars <- c("SNP", "chr", "pos", "pos_b38", 
  "effect_allele", "other_allele", "eaf", 
  "beta", "se", "pval", 
  "Direction", "het_Q", "het_pval", "n_studies", "samplesize")
	  
####################################################################################
##            Extract and combine study-specific sample size info                 ##
####################################################################################

# Function to read study-specific GWAS sum stats
extr_info <- function(index, study, rel_path, prefix, strata = NULL) {
  
  ## Pheno index
  ind <- index
  
  ## Study phenotype
  var_col <- grep(study, names(varlist), value = T) 
  pheno_name <- select(varlist, Variable_name) %>% .[ind,]

  ## Read file and keep only variants with rsIDs
  file <- list.files(path = paste0(sumdat_dir, rel_path), 
	  pattern = prefix,
	  full.names = T) %>%
	  fread %>% 
	  filter(grepl("^rs", SNP)) 
  
  ## Select SNP info
  if(pheno_name %in% bin.out) {
	dt <- select(file, all_of(bin.vars))
  } else if(pheno_name %in% cont.out) {
	dt <- select(file, all_of(cont.vars))
  }
  
  ## Add study to var names
  dt <- rename_at(dt, vars(-SNP), ~ paste0(., '_', study))
  
  print(head(dt))
  return(dt)
 }
  
# Function to read & combine all GWAS sum stats
comb_info <- function(x) {

  ## Pheno
  print(varlist[x, ])
  
  ## Create empty dataframe
  df <- tibble()
  
  ## MOBA
  if(!is.na(varlist$var_MOBA[x])){
    df_moba <- extr_info(
	  index = x,
      study = "MOBA",
      rel_path = "MOBA/mothers/final",
      prefix = paste0("^moba_100k_mum_", varlist$var_MOBA[x])
	) 
	df <- df_moba
  }
  
  ## ALSPAC
  if(!is.na(varlist$var_ALSPAC[x])){
    df_alspac <- extr_info(
	  index = x,
      study = "ALSPAC",
      rel_path = "ALSPAC/mothers/final",
      prefix = paste0("^alspac_mums_", varlist$var_ALSPAC[x])
	) 
	df <- merge(df, df_alspac, by = "SNP", all = T)
  }
  
  ## BiB - white Europeans
  if(!is.na(varlist$var_BIB.WE[x])){
    df_bib.we <- extr_info(
	  index = x,
      study = "BIB.WE",
      rel_path = "BIB-WE/mothers",
      prefix = paste0("^BIB_mother.european_", varlist$var_BIB.WE[x], "_2024"),
	  strata = "WE"
	) 
	df <- merge(df, df_bib.we, by = "SNP", all = T)
  }
  
  ## BiB - south Asians
  if(!is.na(varlist$var_BIB.SA[x])){
    df_bib.sa <- extr_info(
	  index = x,
      study = "BIB.SA",
      rel_path = "BIB-SA/mothers",
      prefix = paste0("^BIB_mother.southasian_", varlist$var_BIB.SA[x], "_2024"),
	  strata = "SA"	
	) 
	df <- merge(df, df_bib.sa, by = "SNP", all = T)
  }  

  ## UKB
  if(!is.na(varlist$var_UKB[x])){
    df_ukb <- extr_info(
	  index = x,
      study = "UKB",
      rel_path = "UKB/mothers/final",
      prefix = paste0("^ukb_mums_", varlist$var_UKB[x], "_2024")
	) 
	df <- merge(df, df_ukb, by = "SNP", all = T)
  }    
  
  ## FinnGen
  if(!is.na(varlist$var_FINNGEN[x])){
    df_finn <- extr_info(
	  index = x,
      study = "FINNGEN",
      rel_path = "FinnGen/mothers",
      prefix = paste0("^FINNGEN-R12[.]", varlist$var_FINNGEN[x])
	) 
	df <- merge(df, df_finn, by = "SNP", all = T)
  }    	
  
  ## Other publicly available GWAS
  if(!is.na(varlist$var_GWAS[x])){
    if(varlist$var_GWAS[x] %in% c("ga_subsamp", "pretb_subsamp", "posttb_subsamp")) {
	   consortium <- "EGG"
	 } else if (varlist$var_GWAS[x] %in% "gdm_subsamp") {
	   consortium <- "GENDIP"
     } else if (varlist$var_GWAS[x] %in% "pe_subsamp") {
	   consortium <- "INTERPREGGEN"
     } else if (varlist$var_GWAS[x] %in% "depr_subsamp") {
	   consortium <- "PGC"
     }
    df_gwas <- extr_info(
	  index = x,
      study = "GWAS",
      rel_path = "Public/mothers",
      prefix = paste0("^", consortium, "_mums_", varlist$var_GWAS[x])
	) 
	df <- merge(df, df_gwas, by = "SNP", all = T)
 }   

  return(df)
}

# Read & combine all GWAS sum stats
n_studies <- map(1:nrow(varlist), ~comb_info(.)) %>%
  set_names(., varlist$Variable_name)

####################################################################################
##                Calculate total sample size for metanalyses                     ##
####################################################################################

# Function to calculate total sample size
sum_n <- function(x) {

  if(varlist$Variable_name[x] %in% bin.out) {
  
	df <- n_studies[[x]] %>%
      mutate(
	    samplesize = rowSums(pick(starts_with("samplesize_")), na.rm = T),
	    ncase = rowSums(pick(starts_with("ncase_")), na.rm = T),
        ncontrol = rowSums(pick(starts_with("ncontrol_")), na.rm = T)
	    ) %>%
	  select(SNP, ncase, ncontrol, samplesize)
	  
  } else if(varlist$Variable_name[x] %in% cont.out) {
  
	df <- n_studies[[x]] %>%
      mutate(
	    samplesize = rowSums(pick(starts_with("samplesize_")), na.rm = T)
	    ) %>%
      select(SNP, samplesize)		
  }
    
  return(df)
} 

## Calculate total sample size, ncase, ncontrol
n_total <- map(1:nrow(varlist), ~sum_n(.)) %>%
  set_names(., varlist$Variable_name)


####################################################################################
##                     Clean and save metanalysis files                           ##
####################################################################################

## Mapping files 

snp_map_dir=paste0(Sys.getenv('MRPREG_sumdat'), "results/liftover/")
snp_map <- fread(paste0(snp_map_dir, "mrpreg_rsids_22m_dbSnp155_hg19_hg38.txt.gz")) %>%
  rename(SNP = rsid, pos = pos_hg19, pos_b38 = pos_hg38) %>%
  select(SNP, chr, pos, pos_b38)
str(snp_map)

## Import METAL file and merge with N information

for (i in 1:nrow(varlist)){
  
  ## Outcome 
    pheno_name <- varlist$Variable_name[i]
    print(pheno_name)
	
  ## Select information file (chr, pos, sample size...)
    n_out <- n_total[[i]]
  
  ## Import METAL file
  	metal_output <- list.files(
  	  path = paste0(sumdat_dir, "Metanalysis/mothers/METAL_original_output/"), 
        pattern = paste0("^metal_", varlist$Variable_name[i], "_1[.]tbl$"),
    	  full.names = T) %>%
  	  fread 
      str(metal_output)
	
  ## Select and rename Variables
   # if(nchar(metal_output$Direction[1])>1) {
    metal_clean <- metal_output %>%
	  mutate(n_studies = HetDf+1) %>%
	  mutate(across(c(Allele1, Allele2), ~toupper(.))) %>%
	  rename(
	    SNP = MarkerName, 
	    effect_allele = Allele1,
		  other_allele = Allele2,
		  eaf = Freq1,
		  beta = Effect,
		  se = StdErr,
	  	pval = 'P-value',
		  het_Q = HetChiSq,
		  het_pval = HetPVal
		) %>%
	  filter(grepl("^rs", SNP)) 
   #  } else {
   # metal_clean <- metal_output %>%
	 # mutate(n_studies = 1) %>%
	 # mutate(across(c(Allele1, Allele2), ~toupper(.))) %>%
	 # rename(
	 #   SNP = MarkerName, 
	 #   effect_allele = Allele1,
	#	  other_allele = Allele2,
	#	  eaf = Freq1,
	#	  beta = Effect,
	#	  se = StdErr,
	 # 	pval = 'P-value'
	#	) %>%
	 # filter(grepl("^rs", SNP))
   #  }
    str(metal_clean)
  
  ## Merge metal and info files
    metal_comb <- merge(metal_clean, n_out, by = "SNP", all.x = TRUE) %>%
	  merge(., snp_map, by = "SNP", all.x = TRUE)
		
  ## Keep if association estimated in > 50% ncase/sample size
	print(paste("N SNPs before QC:", nrow(metal_comb)))
	if(pheno_name %in% bin.out) {
	  metal_comb <- filter(metal_comb, ncase > max(metal_comb$ncase, na.rm = T)/2)
	} else if(pheno_name %in% cont.out) {
	  metal_comb <- filter(metal_comb, samplesize > max(metal_comb$samplesize, na.rm = T)/2)
	}
  	print(paste("N SNPs after removing <= 50% ncase/sample size", nrow(metal_comb)))

  ## Keep if MAF > 0.01
	metal_comb <- filter(metal_comb, eaf > 0.01 & eaf < 0.99)
	print(paste("N SNPs after removing MAF < 1%:", nrow(metal_comb)))

  ## Keep if SE > 0 
	metal_comb <- filter(metal_comb, se > 0)
	print(paste("N SNPs after removing SE = 0:", nrow(metal_comb)))

  ## Keep if heterogenity p-value >  0.001
	#if(max(metal_comb$n_studies)>1) {
   metal_comb <- filter(metal_comb, het_pval > 0.001)
	 print(paste("N SNPs after removing het pval <= 0.001:", nrow(metal_comb)))
  #}
  
  ## Select relevant variables
    if(pheno_name %in% bin.out) {
	  meta.bin.vars <- c(meta.vars, "ncase", "ncontrol")
	  metal_comb <- select(metal_comb, all_of(meta.bin.vars))
    } else if(pheno_name %in% cont.out) {
	  metal_comb <- select(metal_comb, all_of(meta.vars))
    }
    str(metal_comb)
	
  ## Save clean files
    write.table(metal_comb, file = gzfile(
	  paste0(sumdat_dir, "Metanalysis/mothers/final/", "metanalyses-R3.mum.", varlist$Variable_name[i], ".txt.gz")), 
	quote = F, sep = " ", row.names = F)
	
	rm(metal_output, metal_clean, metal_comb)
}
    
q("no")


