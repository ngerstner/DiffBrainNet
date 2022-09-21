##################################################
## Project: DexStim Mouse Brain
## Date: 16.09.2022
## Author: Nathalie
##################################################
# Test if Tcf4 DN is enriched for genes that are actually 
# dysregulated upon Tcf4 mutant mice: 
# (https://pubmed.ncbi.nlm.nih.gov/32015540/)


library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(fgsea)
library(biomaRt)
library(ggplot2)

basedir <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
folder_table <- paste0(basedir,"/tables")
file_tcf4 <- paste0(basedir, "data/reviews/TCF4 PFC DN.xlsx")
file_mutant_tcf4 <- paste0(basedir, "data/reviews/Phan et al_Nat Neuro_2020 TCF4.xls")
file_background <- paste0(folder_table, "/06_background.txt")


## 0. Set up BioMart for ID mapping
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# --> this is a workaround as the normal commands lead to errors


## 1. Read background mouse genes from our dataset
background <- read.table(file_background, 
                         row.names = NULL)


## 2. Read Tcf4 DN from PFC
tcf4_dn <- read_excel(file_tcf4, sheet = "Sheet1")
tcf4_dn <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), # no difference between v3 and v4
                 filters = "mgi_symbol",
                 values = tcf4_dn$`unique identifiers`, mart = mouse)


## 3. Read mutant Tcf4 mice DE genes
mutant <- list()
# adult mice
tcf4_mutant <- read_excel(file_mutant_tcf4, sheet = "Supp Table 2-Adult DESeq2") 
tcf4_demutant <- tcf4_mutant$...1[tcf4_mutant$padj <= 0.05]
tcf4_demutant <- tcf4_demutant[tcf4_demutant %in% background$V1]
mutant[["tcf4_mutant_adult"]] <- tcf4_demutant

# p1 mice
tcf4_mutant <- read_excel(file_mutant_tcf4, sheet = "Supp Table 2-p1 DESeq2") 
tcf4_demutant <- tcf4_mutant$...1[tcf4_mutant$padj <= 0.05]
tcf4_demutant <- tcf4_demutant[tcf4_demutant %in% background$V1]
mutant[["tcf4_mutant_p1"]] <- tcf4_demutant

# both
tcf4_mutant <- read_excel(file_mutant_tcf4, sheet = "Supp Table 2-Both DESeq2")
tcf4_demutant <- tcf4_mutant$...1[tcf4_mutant$padj_p1 <= 0.05 & tcf4_mutant$padj_Adult]
tcf4_demutant <- tcf4_demutant[tcf4_demutant %in% background$V1]
mutant[["both"]] <- tcf4_demutant

## 4. Overrepresentation test
# Check if mutant Tcf4 DE genes are enriched in our network  
foraRes <- fora(pathways = mutant, 
                genes = tcf4_dn$ensembl_gene_id,
                universe = background$V1)


## 5. Write to file
write.csv(as.matrix(foraRes),
          file = paste0(basedir, 
                        "scripts_manuscript/04_PlotsManuscript/Revision/19_Reviewer2_3_Tcf4_enrichments.csv"))

