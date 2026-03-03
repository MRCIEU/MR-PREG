library(data.table)
library(dplyr)
library(purrr)

wdir <- paste0(Sys.getenv('MRPREG_sumdat'), "")

# File info
f_info <- readxl::read_excel(
  path = paste0(wdir, "scripts/GWAS/FinnGen/mothers/", "finngen_outcomes_R12.xlsx"),
  sheet = 1
)

# Mapping file
snp_map_dir=paste0(Sys.getenv('MRPREG_sumdat'), "results/liftover/")
snp_map <- fread(paste0(snp_map_dir, "mrpreg_rsids_22m_dbSnp155_hg19_hg38.txt.gz")) %>%
  rename(SNP = rsid, pos = pos_hg19) %>%
  select(SNP, pos)
str(snp_map)

# Format data function
clean_finn_dat <- function(finn_outname, ncase, ncontrol, mrpreg_outname){
  gc()
  print(finn_outname)
  # List raw file
    f_name <- list.files(
      path = paste0(wdir, "data/GWAS/FinnGen/mothers/R12/"),
	  pattern = paste0("^summary_stats_finngen_R12_",finn_outname),
	  full.name = T
	  )
  # Import raw data
    f <- fread(f_name) 
	print("Raw data:")
	print(head(f))
  # Select relevant columns
    dat <- f %>%
	  rename(
	    SNP = rsids,
		chr = "#chrom",
		pos_b38 = pos,
		effect_allele = alt,
		other_allele = ref,
		eaf = af_alt,
		eaf_cases = af_alt_cases,
		eaf_controls = af_alt_controls,
		se = sebeta
		) %>%
		select(SNP,chr,pos_b38,effect_allele,other_allele,eaf,eaf_cases,eaf_controls,beta,se,pval) %>%
		mutate(samplesize = ncase+ncontrol, ncase = ncase, ncontrol = ncontrol, units = "logodds")
  # Remove MAF < 1%
    dat <- filter(dat, eaf > 0.01 & eaf < 0.99) 
  # Merge summay data with mapping file
    dat <- merge(snp_map, dat, by = "SNP")
  # Convert missing rsIDs in NA
    dat$SNP[dat$SNP==""]<-"NA"
  # Save formatted data  
    print("Formatted data:")
	print(head(dat))
    write.table(dat, file = gzfile(
	  paste0(wdir, "results/GWAS/FinnGen/mothers/", "FINNGEN-R12.",mrpreg_outname,".20241114.txt.gz")
	  ), quote = F, sep = " ", row.names = F)
}
 
# Format data
pwalk(f_info, clean_finn_dat)
#clean_FinnGen("O15_PREG_PROLONGED",6062,200533,"posttb_all")
