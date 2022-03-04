##################################################
## Project: DexStim Mouse Brain
## Date: 16.10.2020
## Author: Nathalie
##################################################
# Parse expression data as input for kimono

library(data.table)
library(dplyr)

regions <- c("AMY", "CER", "PFC", "PVN", "vDG", "dDG", "vCA1", "dCA1")
dex <- c(0,1)

basepath <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
folder_table <- paste0(basepath,"tables/")
output_file <- paste0(basepath,"data/kimono_input/expression_mm_")
files <- list.files(path = folder_table, pattern = "deseq2_expression_vsd.txt$", full.names = TRUE)

# 1. Read expression table for each region
expr_list <- lapply(files, function(x) as.data.frame(t(as.matrix(fread(x),rownames=1))))

# 2. Merge expression tables and remove columns with NAs
expr_all <- bind_rows(expr_list) %>%
  dplyr::select(where(~!any(is.na(.)))) # maybe change this and set NAs to 0

# 3a. Write expression table to file
for (reg in regions){
  for (d in dex){
    dex_str <- ifelse(d == 0, "CNTRL", "DEX")
    phen <- expr_all[startsWith(row.names(expr_all),reg),]
    phen <- phen[grepl(dex_str, row.names(phen)),]
    fwrite(phen, file = paste0(output_file, reg, "_dex", d, ".csv"), row.names = TRUE)
  }
}

# 3b. Scale all expression values together and write to file
expr_scaled <- as.data.frame(scale(expr_all))
# expr_scaled[,1:10]
# expr_all[,1:10]
for (reg in regions){
  for (d in dex){
    dex_str <- ifelse(d == 0, "CNTRL", "DEX")
    phen <- expr_scaled[startsWith(row.names(expr_scaled),reg),]
    phen <- phen[grepl(dex_str, row.names(phen)),]
    fwrite(phen, file = paste0(output_file, reg, "_dex", d, "_scaled.csv"), row.names = TRUE)
  }
}
