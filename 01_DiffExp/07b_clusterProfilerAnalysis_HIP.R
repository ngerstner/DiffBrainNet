##################################################
## Project: DexStim Mouse Brain
## Date: 19.04.2021
## Author: Nathalie
##################################################
# Functional annotation for HIP with clusterProfiler
# make figure for manuscript

library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(data.table)


basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"

regions <-
  c("vDG", "dDG", "vCA1", "dCA1")



# 1. Read DE tables from HIP regions ----------

list_reg_sig <- list()
background <- NULL

for (reg in regions) {
  res <-
    fread(
      file = paste0(
        basepath,
        "tables/02_",
        reg,
        "_deseq2_Dex_1_vs_0_lfcShrink.txt"
      ),
      sep = "\t"
    )
  na_indices <- which(is.na(res$padj))
  res$padj[na_indices] <- 1
  res_sig <- res[res$padj <= 0.1, ]
  # res_sig <- res[res$log2FoldChange >= 1]
  list_reg_sig[[reg]] <- res_sig
  background <- res$Ensembl_ID
}



# 2. Concatenate DE tables -----------------

data <- bind_rows(list_reg_sig, .id = "region")



# 3. GO enrichment -------------------------

# IMPORTANT: which background?

for (reg in regions){
  
  genes <- list_reg_sig[[reg]]$Ensembl_ID
  # background <- unique(data$Ensembl_ID)

  # TODO: decide on maxGSSize --> with 10000 very similar results to anRichment
  ego <- enrichGO(gene          = genes,
                  universe      = background,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  minGSSize = 10,    # min number of genes associated with GO term
                  maxGSSize = 10000, # max number of genes associated with GO term
                  readable      = TRUE)
  print(head(ego, n = 20))

}
