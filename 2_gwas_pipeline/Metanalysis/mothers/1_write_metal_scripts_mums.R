# Packages
library(dplyr)

# Directories
w_dir=paste0(Sys.getenv('MRPREG_sumdat'))
scripts_dir=dat_dir=paste0(w_dir, "/scripts/GWAS/Metanalysis/mothers/")
code_dir="./Linux/"
dat_dir=paste0(w_dir, "/results/GWAS/")
out_dir=paste0(dat_dir, "Metanalysis/mothers/")

# List of MR-PREG outcomes
varlist <- read.csv(paste0(scripts_dir, "vars_metanalysis_v5_mum.csv"), header = T)

# Sample size for each GWAS
n_gwas <- read.table(paste0(dat_dir, "samplesize/N_gwas.txt"), header=T)
low_ncase <- filter(n_gwas, ncase < 50) %>%  
  distinct(study, outcome)
  
# Remove study outcomes that should not contribute to metanalyses
lowncase_studies <- unique(low_ncase$study)
for(s in lowncase_studies) {
  print(s)
  # Outcomes with ncase < 50 in the study (mums or children)
  excl_out <- filter(low_ncase, study == s) %>% pull(outcome)
  print(excl_out)
  # Column name for study
  if(grepl("BIB", s)) {
   s <- sub("-", ".", s)
  }
  var_col <- paste0("var_", s)
  # Create logical vector where outcome to be excluded == TRUE
  to_excl <- varlist[,var_col] %in% excl_out
  # Turn outcome to be excluded into NA
  varlist[,var_col] <- ifelse(to_excl, NA, varlist[,var_col])
}
write.table(varlist, file = paste0(out_dir, "vars_metanalysis_v5_mum_final.txt"), quote = F, row.names = F, sep="\t")

# Write METAL scripts for each outcome
for(i in 1:nrow(varlist)){
  
  ## Print variable info
  print(varlist[i,])
  
  ## Metal script file
  metal_script=paste0(code_dir, "metal_", varlist$Variable_name[i], ".txt")
  
  ## Header
  sink(metal_script)
  cat("GENOMICCONTROL ON")
  cat("\n")
  cat("AVERAGEFREQ ON")
  cat("\n")
  cat("SEPARATOR WHITESPACE")
  cat("\n")
  cat("SCHEME STDERR")
  cat("\n")
  cat("MARKER SNP")
  cat("\n")
  cat("ALLELE effect_allele other_allele")
  cat("\n")
  cat("FREQ eaf")
  cat("\n")
  cat("EFFECT beta")
  cat("\n")
  cat("STDERR se")
  cat("\n")
  sink()
  
  ## MoBa
  if(!is.na(varlist$var_MOBA[i])){
   ### Find input file
    file_moba <- list.files(path = paste0(dat_dir, "MOBA/mothers/final"), 
	  pattern = paste0("^moba_100k_mum_", varlist$var_MOBA[i]),
	  full.names = T
	)
	print(file_moba)
   ### Process input file
    sink(metal_script, append = T)
    cat(paste0("PROCESS ", file_moba))
    cat("\n")
    sink()
  }
  
  ## ALSPAC
  if(!is.na(varlist$var_ALSPAC[i])){
   ### Find input file
    file_alspac <- list.files(path = paste0(dat_dir, "ALSPAC/mothers/final"), 
	  pattern = paste0("^alspac_mums_", varlist$var_ALSPAC[i]),
	  full.names = T
	  )
	print(file_alspac)
   ### Process input file
	sink(metal_script, append = T)
	cat(paste0("PROCESS ", file_alspac))
	cat("\n")
	sink()
  }
  
  ## BiB - white Europeans
  if(!is.na(varlist$var_BIB.WE[i])){
   ### Find input file
    file_bib_we <- list.files(path = paste0(dat_dir, "BIB-WE/mothers"), 
      pattern = paste0("^BIB_mother.european_", varlist$var_BIB.WE[i], "_2024"),
	  full.names = T
	)
	print(file_bib_we)
   ### Process input file
    sink(metal_script, append = T)
	cat(paste0("PROCESS ", file_bib_we))
	cat("\n")
	sink()
  }  
  
  ## BiB - south Asians
  if(!is.na(varlist$var_BIB.SA[i])){
   ### Find input file
    file_bib_sa <- list.files(path = paste0(dat_dir, "BIB-SA/mothers"), 
      pattern = paste0("^BIB_mother.southasian_", varlist$var_BIB.SA[i], "_2024"),
	  full.names = T
	)
	print(file_bib_sa)
   ### Process input file
	sink(metal_script, append = T)
	cat(paste0("PROCESS ", file_bib_sa))
	cat("\n")
	sink()	
  }
   
  ## FinnGen
  if(!is.na(varlist$var_FINNGEN[i])){
   ### Find input file
    file_finngen <- list.files(path = paste0(dat_dir, "FinnGen/mothers"), 
      pattern = paste0("^FINNGEN-R12[.]", varlist$var_FINNGEN[i]),
	  full.names = T
	)	
	print(file_finngen)
   ### Process input file
    sink(metal_script, append = T)
	cat(paste0("PROCESS ", file_finngen))
	cat("\n")
	sink()	
  }
  
  ## UKB
  if(!is.na(varlist$var_UKB[i])){
   ### Find input file
    file_ukb <- list.files(path = paste0(dat_dir, "UKB/mothers/final"), 
      pattern = paste0("^ukb_mums_", varlist$var_UKB[i], "_2024"),
	  full.names = T
	)	
	print(file_ukb)
	### Process input file
    sink(metal_script, append = T)
	cat(paste0("PROCESS ", file_ukb))
	cat("\n")
	sink()	
  }
  
  ## Publicly available GWAS
  if(!is.na(varlist$var_GWAS[i])){
    ### Find input file
	 if(varlist$var_GWAS[i] %in% c("ga_subsamp", "pretb_subsamp", "posttb_subsamp")) {
	   consortium <- "EGG"
	   } else if (varlist$var_GWAS[i] %in% "gdm_subsamp") {
	   consortium <- "GENDIP"
     } else if (varlist$var_GWAS[i] %in% "pe_subsamp") {
	   consortium <- "INTERPREGGEN"
     } else if (varlist$var_GWAS[i] %in% "depr_subsamp") {
	   consortium <- "PGC"
     }
     file_pub <- list.files(path = paste0(dat_dir, "Public/mothers"), 
       pattern = paste0("^", consortium, "_mums_", varlist$var_GWAS[i]),
	   full.names = T
	 )
	 print(file_pub)
    ### Process input file
    sink(metal_script, append = T)
    cat(paste0("PROCESS ", file_pub))
	cat("\n")
    sink()
  }
  
  ## Tail
  sink(metal_script, append = T)
  cat(paste0("OUTFILE ", out_dir, "METAL_original_output/metal_", varlist$Variable_name[i], "_ .tbl"))
  cat("\n")
  cat("ANALYZE HETEROGENEITY")
  cat("\n")
  cat("QUIT")
  sink()
}

q("no")
