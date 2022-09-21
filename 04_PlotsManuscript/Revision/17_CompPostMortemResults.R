##################################################
## Project: DexStim Mouse Brain
## Date: 14.09.2022
## Author: Nathalie
##################################################
# Compare DEs and DNs with DE genes from post mortem study
# (https://pubmed.ncbi.nlm.nih.gov/30545856/)


library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(biomaRt)
library(ggplot2)
library(UpSetR)
library(rlist)


# define pathes
basedir <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse"
folder_table <- paste0(basedir,"/tables")
file_pm <- paste0(basedir, "/data/reviews/PostmortemBrain/aat8127_table_s1.xlsx")
files_de <- list.files(path = folder_table, pattern = "deseq2_Dex_0_vs_1_lfcShrink.txt$", full.names = TRUE)
print(sub(".*02_(\\w*)_deseq2.*","\\1",files_de))
files_dn <- list.files(path = paste0(folder_table, "/coExpression_kimono/03_AnalysisFuncoup/"), 
                       pattern = "_funcoup_treatment_nodebetweennessNorm_betacutoff0.01.csv$", full.names = TRUE)
regions_files <- sub(".*03a_(\\w*)_funcoup_treatment_nodebetween.*","\\1",files_dn)
file_background <- paste0(folder_table, "/06_background.txt")

## 1. Set up BioMart
# map mouse Ensembl Ids to human Ensembl Ids --> homology mapping
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# --> this is a workaround as the normal commands lead to errors


## 2. Read background genes and get human background IDs
background <- read.table(file_background, 
                         row.names = NULL)
background_human <- getLDS(attributes=c("ensembl_gene_id"),
                           filters="ensembl_gene_id", values=background$V1, mart=mouse,
                           attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1


## 3. Read differentially expressed genes and transcripts from postmortem study
list_pm_genes <- list()
for (r in c("DGE", "DTE", "DTU")){
  # read file including ASD, SCZ and BD results
  data <- read_excel(file_pm, sheet=r) 
  
  list_pm_genes[[r]] <- list()
  
  for (d in c("ASD", "BD", "SCZ")){
    data_d <- data %>%
      rename_with(~ gsub('DTU\\.', '', .x)) %>% # DTU list has prefix "DTU" in column names
      dplyr::select("ensembl_gene_id", starts_with(d)) %>%
      rename_with(~ gsub(paste0(d,'\\.'), '', .x)) %>%
      rename_with(tolower) %>% # DTU table has FDR instead of fdr
      dplyr::filter(fdr <= 0.05)
    list_pm_genes[[r]][[d]] <- data_d$ensembl_gene_id[data_d$ensembl_gene_id %in% background_human]
    }
}


## 4. Read DE tables from all regions 
expr_list <- lapply(files_de, function(x) fread(x) %>% filter(padj <= 0.1))
names(expr_list) <- regions_files

de_list <- lapply(expr_list, function(x) x$V1)
de_list_human <- lapply(de_list, 
                        function(x) getLDS(attributes=c("ensembl_gene_id"),
                                           filters="ensembl_gene_id", values=x, mart=mouse,
                                           attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1)


## 5. Upset plot of PFC and post mortem brain DE genes
list_flat <- list("PFC.mouse" = de_list_human[["PFC"]],
               "PFC.human.ASD" = list_pm_genes[["DGE"]][["ASD"]],
               "PFC.human.BD" = list_pm_genes[["DGE"]][["BD"]],
               "PFC.human.SCZ" = list_pm_genes[["DGE"]][["SCZ"]])
pdf(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/17_Reviewer2_2_DEGpostmortem.pdf"),
    width = 12, height = 8)
par(mar=c(6,5,4,2) + 0.1)
print(upset(fromList(list_flat), nsets = 4, order.by = "freq",
            text.scale = c(1.8, 1.9, 1.8, 1.9, 1.9, 1.9),
            sets.x.label = "#DE genes",
            mainbar.y.label = "#DE genes in intersection"))
dev.off()
