##################################################
## Project: DexStim Mouse Brain
## Date: 16.10.2020
## Author: Nathalie
##################################################
# Mapping between genes and phenotypes as input for kimono
# Sufficient to do this once, not for each region and dex status

library(data.table)
library(dplyr)
library(tidyr)

basepath <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
folder_table <- paste0(basepath,"tables/")
output_file <- paste0(basepath,"data/kimono_input/prior_expr_bio_mm_regDex.csv")
pheno_file <- paste0(basepath,"data/kimono_input/phenotypes_mm_AMY_dex0.csv")
expr_file <- paste0(basepath,"data/kimono_input/expression_mm_AMY_dex0.csv")


# 1. Read phenotype and expression file
pheno <- fread(pheno_file)
expr <- fread(expr_file)


# 2. Create all pairs of pheno and gene
map <- tidyr::crossing("gene" = colnames(expr)[-1], "bio" = colnames(pheno)[-1])


# 3. Write mapping to file 
fwrite(map, output_file)
