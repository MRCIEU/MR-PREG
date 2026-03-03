## Script to calculate PCs 


 

## Pkgs
#remotes::install_github("danjlawson/pcapred")
#remotes::install_github("danjlawson/pcapred.largedata")
library(pcapred)
library(pcapred.largedata)
#install.packages('R.utils')
#library(R.utils)
library(dplyr)

setwd("")

## Ref
ref=readreference(ukb_pcs_200())

 
#dat=readbed("merged-pcs") 































## CE - SA

dat=readbed("merged-pcs-ce-sa") 

# Find the index of the 'snpstats' data frame within the 'ref' list
snpstats_index <- which(names(ref) == "snpstats")

# Check if 'snpstats' was found in the list
if (length(snpstats_index) == 0) {
  cat("The 'snpstats' data frame was not found in the 'ref' list.\n")
} else {
  # Access the 'snpstats' data frame and add the 'markerID' column by pasting info
  snpstats <- ref[[snpstats_index]]
  
  # Add the 'mID' column using the mutate() function
  snpstats <- snpstats %>%
    mutate(ID = paste(X.CHROM, position, REF, ALT, sep = ":"))
  
  # Assign the modified data frame back to the list
  ref[[snpstats_index]] <- snpstats
}

rm(snpstats)

## Merge
dat=mergeref(dat,ref=ref)
dat$mode<-"flashpca"
## Runs to here
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
writepred("saige-pcs-pcapred-eigenvals-ce-sa",dat,pred) 

## GSA - WE

dat=readbed("merged-pcs-gsa-we") 

# Find the index of the 'snpstats' data frame within the 'ref' list
snpstats_index <- which(names(ref) == "snpstats")

# Check if 'snpstats' was found in the list
if (length(snpstats_index) == 0) {
  cat("The 'snpstats' data frame was not found in the 'ref' list.\n")
} else {
  # Access the 'snpstats' data frame and add the 'markerID' column by pasting info
  snpstats <- ref[[snpstats_index]]
  
  # Add the 'mID' column using the mutate() function
  snpstats <- snpstats %>%
    mutate(ID = paste(X.CHROM, position, REF, ALT, sep = ":"))
  
  # Assign the modified data frame back to the list
  ref[[snpstats_index]] <- snpstats
}

rm(snpstats)

## Merge
dat=mergeref(dat,ref=ref)
dat$mode<-"flashpca"
## Runs to here 
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
writepred("saige-pcs-pcapred-eigenvals-gsa-we",dat,pred)

## GSA - SA 

dat=readbed("merged-pcs-gsa-sa") 

# Find the index of the 'snpstats' data frame within the 'ref' list
snpstats_index <- which(names(ref) == "snpstats")

# Check if 'snpstats' was found in the list
if (length(snpstats_index) == 0) {
  cat("The 'snpstats' data frame was not found in the 'ref' list.\n")
} else {
  # Access the 'snpstats' data frame and add the 'markerID' column by pasting info
  snpstats <- ref[[snpstats_index]]
  
  # Add the 'mID' column using the mutate() function
  snpstats <- snpstats %>%
    mutate(ID = paste(X.CHROM, position, REF, ALT, sep = ":"))
  
  # Assign the modified data frame back to the list
  ref[[snpstats_index]] <- snpstats
}

rm(snpstats)

## Merge
dat=mergeref(dat,ref=ref)
dat$mode<-"flashpca"
## Runs to here
pred=predictpcs(dat)  # Predict the first 18 UK Biobank PCs
writepred("saige-pcs-pcapred-eigenvals-gsa-sa",dat,pred)