## Merging both sets PCs to each other and outcome data ##






## Pull file paths from config.sh

config_lines <- readLines("config.sh")















dir<-grep("GWAS_dat=", config_lines, value = TRUE)
dir<-gsub("GWAS_dat=", "", dir)
dir<-gsub("^\\s+|\\s+$|\"|'", "", dir)

library(dplyr)
library(data.table)

## Outcome data - basic qc - rename ID column and add extra IDD/FID so SAIGE runs


outcomes <- read.table(paste0(dir,"/bib_outcomes_master_110624.txt"),header=T)
outcomes <- outcomes %>% rename(IID = lab_id)
outcomes <- outcomes %>% mutate(FID = IID)

## Read PCs and add IID column


externals.we <- read.table(paste0(dir,"/white_mothers_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_projections_20.txt"),header=T)
externals.sa <- read.table(paste0(dir,"/pakistani_mothers_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_projections_20.txt"),header=T)
externals.we <- externals.we %>% mutate(IID = FID)
externals.sa <- externals.sa %>% mutate(IID = FID)

head(externals.we)
head(externals.sa)

## Read other PCs sets 



internals.we.ce <- read.table(paste0(dir,"/saige-pcs-pcapred-eigenvals"))
internals.we.gsa <- read.table(paste0(dir,"/saige-pcs-pcapred-eigenvals-gsa-we"))
internals.sa.ce <- read.table(paste0(dir,"/saige-pcs-pcapred-eigenvals-ce-sa"))
internals.sa.gsa <- read.table(paste0(dir,"/saige-pcs-pcapred-eigenvals-gsa-sa"))

## Make adjustments - the col names are strange - remove # and change 
## PC names to PC


# Function to apply your code and select the first 22 columns
process_dataframe <- function(df) {
  df %>%
    select(-1) %>%
    rename(FID = V2) %>%
    rename_at(vars(starts_with("V")), ~ paste0("PC", seq_len(length(.)))) %>%
    select(1:21) %>%
	mutate(IID = FID) %>%
	select(1, IID, everything())
}

# Apply the function to all dataframes and save the results
internals.we.ce <- process_dataframe(internals.we.ce)
internals.we.gsa <- process_dataframe(internals.we.gsa)
internals.sa.ce <- process_dataframe(internals.sa.ce)
internals.sa.gsa <- process_dataframe(internals.sa.gsa)

## Check overlap
a <- intersect(internals.we.ce$FID, externals.we$FID)
b <- intersect(internals.sa.ce$FID, externals.sa$FID)

## All fine 

c <- intersect(internals.we.gsa$FID, externals.we$FID)
d <- intersect(internals.sa.gsa$IID, externals.sa$IID)

## Merge chips 

internals.we <- rbind(internals.we.ce,internals.we.gsa)
internals.sa <- rbind(internals.sa.ce,internals.sa.gsa)

## Merge PCs


white.pcs <- inner_join(internals.we,externals.we,by = c("IID","FID"))
sa.pcs <- inner_join(internals.sa,externals.sa,by = c("IID","FID"))

# Get a list of all data frames in the workspace
df_list <- mget(ls(pattern = "pcs"))

# Save dataframes in df_list
for (name in names(df_list)) {
  file_name <- paste0(name, ".master.txt")
  write.table(df_list[[name]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}

## merge w pheno 

pheno.we <- inner_join(outcomes,white.pcs,by = "IID")
pheno.sa <- inner_join(outcomes,sa.pcs,by = "IID")

# Get a list of all data frames in the workspace
df_list <- mget(ls(pattern = "^pheno\\."))

# Save dataframes in df_list
for (name in names(df_list)) {
  file_name <- paste0(name, ".master.txt")
  write.table(df_list[[name]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}

#######################################################################################



outcomes <- read.table(paste0(dir,"bib_outcomes_master_child_110624.txt"),header=T)
outcomes <- outcomes %>% rename(IID = lab_id)
outcomes <- outcomes %>% mutate(FID = IID)

## Read PCs and add IID column

externals.we <- read.table(paste0(dir,"white_children_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_projections_20.txt"),header=T)
externals.sa <- read.table(paste0(dir,"pakistani_children_unrelated_maf_0_01_geno_0_03_no_long_range_ld_pruned_projections_20.txt"),header=T)
externals.we <- externals.we %>% mutate(IID = FID)
externals.sa <- externals.sa %>% mutate(IID = FID)

head(externals.we)
head(externals.sa)

internals.we.ce <- read.table(paste0(dir,"saige-pcs-pcapred-eigenvals"))
internals.we.gsa <- read.table(paste0(dir,"saige-pcs-pcapred-eigenvals-gsa-we"))
internals.sa.ce <- read.table(paste0(dir,"saige-pcs-pcapred-eigenvals-ce-sa"))
internals.sa.gsa <- read.table(paste0(dir,"saige-pcs-pcapred-eigenvals-gsa-sa"))

# Function to apply your code and select the first 22 columns
process_dataframe <- function(df) {
  df %>%
    select(-1) %>%
    rename(FID = V2) %>%
    rename_at(vars(starts_with("V")), ~ paste0("PC", seq_len(length(.)))) %>%
    select(1:21) %>%
	mutate(IID = FID) %>%
	select(1, IID, everything())
}

# Apply the function to all dataframes and save the results
internals.we.ce <- process_dataframe(internals.we.ce)
internals.we.gsa <- process_dataframe(internals.we.gsa)
internals.sa.ce <- process_dataframe(internals.sa.ce)
internals.sa.gsa <- process_dataframe(internals.sa.gsa)

## Check overlap
a <- intersect(internals.we.ce$FID, externals.we$FID)
b <- intersect(internals.sa.ce$FID, externals.sa$FID)
c <- intersect(internals.we.gsa$FID, externals.we$FID)
d <- intersect(internals.sa.gsa$IID, externals.sa$IID)

## Combine chips by ancetsry 
internals.we <- rbind(internals.we.ce,internals.we.gsa)
internals.sa <- rbind(internals.sa.ce,internals.sa.gsa)

## Combine pc's by ancestry
white.pcs <- inner_join(internals.we,externals.we,by = c("IID","FID"))
sa.pcs <- inner_join(internals.sa,externals.sa,by = c("IID","FID"))

# Get a list of all data frames in the workspace
df_list <- mget(ls(pattern = "pcs"))

# Save dataframes in df_list
for (name in names(df_list)) {
  file_name <- paste0(name, ".children.txt")
  write.table(df_list[[name]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}

## Combine w pheno 

pheno.we <- inner_join(outcomes,white.pcs,by = "IID")
pheno.sa <- inner_join(outcomes,sa.pcs,by = "IID")

## Save

# Get a list of all data frames in the workspace
df_list <- mget(ls(pattern = "^pheno\\."))

# Save dataframes in df_list
for (name in names(df_list)) {
  file_name <- paste0(name, ".master.children.txt")
  write.table(df_list[[name]], file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}





