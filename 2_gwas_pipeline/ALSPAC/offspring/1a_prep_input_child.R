############################################################################
#                                Set - UP                                  #
############################################################################

##
## Remove (to clear the work environment)
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
## Check what packages are out of date
#old.packages()

##
## Updating packages (update packages already installed)
#update.packages(ask = "FALSE")

##
## Install packages
#install.packages(c("data.table", "purrr", dplyr"))

##
## Load library
library(dplyr)
library(purrr)
library(data.table)

##
## Data directory
gen.dir <- paste0(Sys.getenv('ALSPAC_HRC'), "/bgen/")
phen.dir <- paste0(Sys.getenv('MRPREG_sumdat'), "data/GWAS/ALSPAC/offspring/")
pcs.dir  <- paste0(Sys.getenv('ALSPAC_GENO'), "/derived/principal_components/plink/children/")
out.dir <- paste0(Sys.getenv('MRPREG_sumdat'), "data/GWAS/ALSPAC/offspring/")


############################################################################
#                 Merge phenotypic dataset to .sample file                 #
############################################################################

##
## Import sample file

samp <- fread(paste0(gen.dir, "data.sample")) %>% select(ID_1:missing, sex)


##
## Import file containing MR-PREG outcomes data

phen <- read.table(paste0(phen.dir, "pheno_mr-preg_2024_11_13.txt")) %>%
  mutate(aln_child = paste0(aln, qlet)) %>%
  select(-aln, -qlet, -aln1, -aln2) 
# P.S.: pheno file created by Ana in RDSF ieu2\p6\112


##
## Create covariables file

covar <- read.table(paste0(pcs.dir, "data.eigenvec")) 
covar_names <- c("FID", "IID", paste0("PC", 1:20))
names(covar) <- covar_names
covar <- covar %>% select(FID:PC10)

		
## 
## Groups of variables

# ID variables
id.var <- grep("^aln", names(phen), value = T)

# Continuous variables
cont.var <- grep(("^zbw|^ga_|^apgar|^bf_dur|^nvp_sev"), names(phen), value = T)

# Variable not derived for GWAS
excl.var <- c("misc_index_all", "misc_i_subsamp", "sb_index_all", "sb_i_subsamp", "bw")

# Binary variables
bin.var <- setdiff(names(phen), c(id.var, cont.var, excl.var))

# All outcomes
outlist <- c(bin.var, cont.var)


##
## Recode outcomes for Plink

recode_func <- function(x) {
  case_match(x, "0" ~ 1, "1" ~ 2)

}

plink_phen <- phen %>%
  mutate_at(bin.var, ~as.character(.)) %>%
  mutate_at(bin.var, ~recode_func(.)) %>%
  rename(IID = aln_child) %>%
  mutate(FID = IID)  %>%
  relocate(FID, IID)
  

##
## Save files

# Outcomes
write.table(plink_phen, file = paste0(out.dir, "alspac_child.phen"), quote = F, sep = " ", row.names = F)

# List of outcomes
write.table(outlist, file = paste0(out.dir, "outlist.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# List of binary outcomes
write.table(bin.var, file = paste0(out.dir, "bin_pheno.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# List of continuous outcomes
write.table(cont.var, file = paste0(out.dir, "cont_pheno.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# Sample files
write.table(samp, file = paste0(out.dir, "alspac_child.sample"), quote = F, row.names = F)	

# Covariables
write.table(covar, file = paste0(out.dir, "alspac_child.pcs"), quote = F, row.names = F)

q("no")
