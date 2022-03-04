##################################################
## Project: DexStim Mouse Brain
## Date: 16.10.2020
## Author: Nathalie
##################################################
# Parse phenotype data as input for kimono

library(data.table)
library(dplyr)
library(tidyr)

basepath <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
folder_table <- paste0(basepath,"tables/")
output_file <- paste0(basepath,"data/kimono_input/phenotypes_mm_")
files <- list.files(path = folder_table, pattern = "_deseq2_bio_variables.txt$", full.names = TRUE)


# 1. Read covariable table for each region
biol_list <- lapply(files, fread)
biol_all <- bind_rows(biol_list)


# 2. Fill empty SV cells with 0
biol_all <- biol_all %>% mutate_all(~replace(., is.na(.), 0))


# 3. Single regions and dex status
regions <- unique(biol_all$Region)
dex <- c(0,1)

# 4. Print table of each region to file
for (reg in regions){
  for (d in dex){
    phen <- biol_all %>% filter(Region == reg, Dex == d) %>%
      dplyr::select(-Region, -Dex)
    fwrite(phen, file = paste0(output_file, reg, "_dex", d, ".csv"))
  }
}
