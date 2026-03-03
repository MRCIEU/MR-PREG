##############################################################################################################################################################
######### script to create separate UKB phenotype files for ###################
######### each phenotype to run UKB GWAS using REGENIE ########################
###############################################################################








### set up ###################################################################
library(dplyr)
rm(list = ls())






scratch <- Sys.getenv("MRPREG_sumdat")
out_dir <- paste0(scratch, "/data/GWAS/UKB/mothers/")
setwd(out_dir)

### read in phenotype file ###################################################
pheno_questionnaire <- read.table("Phenotype_MRPREG.txt", header = TRUE)

pheno_mat_records <- read.delim("phenotypic_data_maternity.txt")

# combine using FID, ID
all_pheno <- left_join(pheno_questionnaire, pheno_mat_records, by = join_by(FID, IID))
dim(pheno_questionnaire); dim(pheno_mat_records); dim(all_pheno)
# leaves maternity record outcomes as NA where these are missing

### create new individual phenotype files #####################################

continuous_variables <- c(
    # from questionnaire
    "zbw_all",
    # from maternity records
    "zbw_all_obs", "zbw_subsamp_obs",
    "ga_all", "ga_subsamp"
    )

binary_variables <- c(
    # from questionnaire
    "lbw_all", "hbw_all", "gdm_all", "gdm_subsamp", "gdm_sr_all",
    "gdm_sr_subsamp", "gdm_hes_all", "gdm_hes_subsamp",
    "hyperemesis_all",
    "hdp_all", "hdp_subsamp", "gh_all", "gh_subsamp", "pe_all",
    "pe_subsamp", "depr_all", "depr_subsamp", "pregloss", "sb_all",
    "sb_subsamp", "misc_all", "misc_subsamp", "s_misc_all",
    "s_misc_subsamp", "r_misc_all", "r_misc_subsamp",
    # from maternity records
    "sga", "lga", "pretb_subsamp", "pretb_all", "vpretb_subsamp", "vpretb_all",
    "posttb_subsamp", "posttb_all", "cs", "em_cs", "el_cs", "induction",
    "lbw_all_obs", "lbw_subsamp_obs", "hbw_all_obs", "hbw_subsamp_obs"
    )
# one too many binary variables compared to spreadsheet -
# added "induction", available in maternity records
# also only "sga" and "lga" in pheno files, not _all and _subsamp

continuous_pheno <- all_pheno[, c("FID", "IID", continuous_variables)]
binary_pheno <- all_pheno[, c("FID", "IID", binary_variables)]

### change name of hyperemesis_all to hyp_all
names(binary_pheno)[names(binary_pheno) == "hyperemesis_all"] <- "hyp_all"
binary_variables <- replace(binary_variables,
    binary_variables == "hyperemesis_all", "hyp_all")

### write out individual phenotype files ######################################

# continuous loop for all
for(i in continuous_variables){
    tmp_cols <- c("FID", "IID", i)
    print(tmp_cols) # to check progress as it runs
    tmp_df <- select(continuous_pheno, all_of(tmp_cols))
    write.table(
        tmp_df,
        paste0(out_dir, "pheno_files/cont_", i, ".txt"),
        quote = FALSE, row.names = FALSE)
}

# binary loop for all
for(i in binary_variables){
    tmp_cols <- c("FID", "IID", i)
    print(tmp_cols) # to check progress as it runs
    tmp_df <- select(binary_pheno, all_of(tmp_cols))
    write.table(
        tmp_df,
        paste0(out_dir, "pheno_files/bin_", i, ".txt"),
        quote = FALSE, row.names = FALSE)
}

### write out .txt lists of binary & continuous phenotypes for REGENIE ########
# table for REGENIE to use
write.table(continuous_variables,
    paste0(out_dir, "cont_pheno.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)

# table for REGENIE to use
write.table(binary_variables,
    paste0(out_dir, "bin_pheno.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)

# table with all for formatting
write.table(c(binary_variables, continuous_variables),
    paste0(out_dir, "all_pheno.txt"),
    quote = FALSE, row.names = FALSE, col.names = FALSE)

q("no")