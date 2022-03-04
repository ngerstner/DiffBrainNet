##################################################
## Project: DexStim Mouse Brain
## Date: 09.10.2020
## Author: Nathalie
##################################################
# Parse expression data to mean value per Region/Dex status

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(org.Mm.eg.db)

folder_table <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables"
output_file <- file.path(folder_table, "12_meanExp_regionDex.csv")
files <- list.files(path = folder_table, pattern = "deseq2_expression_vsd.txt$", full.names = TRUE)

# 1. Read expression table for each region
expr_list <- lapply(files, function(x) as.data.frame(t(as.matrix(fread(x),rownames=1))))

# 2. Merge expression tables and remove columns with NAs
expr_all <- bind_rows(expr_list) %>%
  dplyr::select(where(~!any(is.na(.)))) %>%   # maybe change this and set NAs to 0
  tibble::rownames_to_column("sample") %>%    # copy rownames to column
  separate(sample, c("region", "dex"), "_")   # separate former row name into region and dex status
expr_all$dex <- str_remove(expr_all$dex, "\\d+")  # remove mouse number of dex status

# 3. Get mean expression per region/dex group
expr_mean <- expr_all %>%
  group_by(region, dex) %>%
  summarise(across(everything(), mean)) %>%   # mean per group for all genes
  mutate(x = paste(region, dex, sep="_")) %>% # concatenate region and dex
  tibble::column_to_rownames("x") %>%         # and use as rownames
  dplyr::select(-region,-dex)                 # remove region and dex column
expr_mean <- as.data.frame(t(expr_mean))

# 4. Add gene symols
expr_mean$gene_symbol <- mapIds(org.Mm.eg.db, keys = rownames(expr_mean), 
                                keytype = "ENSEMBL", column="SYMBOL")

# 4. Write mean values to file
fwrite(expr_mean, file = output_file, row.names = TRUE)
