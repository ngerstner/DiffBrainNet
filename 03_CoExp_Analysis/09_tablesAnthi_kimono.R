##################################################
## Project: DexStim Mouse Brain
## Date: 15.02.2020
## Author: Nathalie
##################################################
# Make tables for Anthi with kimono results

setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

library(org.Mm.eg.db)
library(dplyr)
library(data.table)
library(anRichment)
library(anRichmentMethods)

regions <- c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1")
mode <- "differential"

# 1. read hub gene tables from all regions ----------

list_reg <- list()
for (reg in regions){
  if (mode == "differential"){
    res <- fread(file=paste0("tables/coExpression_kimono/03_AnalysisFuncoup/", 
                             "04_", reg, "_funcoup_", mode, "_nodebetweennessNorm_betacutoff0.01.csv"))
  } else {
    res <- fread(file=paste0("tables/coExpression_kimono/03_AnalysisFuncoup/", 
                             "03_", reg, "_funcoup_", mode, "_nodebetweennessNorm_betacutoff0.01.csv"))
  }
  res <- res %>%
    filter(nodebetweenness_norm >= 1)
  res$entrez <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                       keytype = "ENSEMBL", column="ENTREZID")
  list_reg[[reg]] <- res
}


# 2. check uniqueness of hub genes ---------------

for (reg in regions){
  index_reg <- which(regions == reg)
  df <- bind_rows(list_reg[-index_reg], .id="region")
  # find regions where gene is also hub
  list_reg[[reg]]$regions_hub <- sapply(list_reg[[reg]]$ensembl_id, 
                                       function(x) paste(df[df$ensembl_id == x,]$region, collapse = " "))
  # boolean if gene is hub uniquely in this region
  list_reg[[reg]]$unique_hub <- sapply(list_reg[[reg]]$regions_hub, 
                                      function(x) x == "")
}


# 3. GO enrichment for the genes of each region ------------------

go_enrichment_all <- function(df_reg, GOcoll, unique, background){
  if (unique){
    genes <- df_reg$entrez[df_reg$unique_hub]
  } else {
    genes <- df_reg$entrez
  }
  modules <- rep("not_significant", length(background))
  modules[which(background %in% genes)] <- "significant"
  
  # enrichment
  GOenrichment <- enrichmentAnalysis(
    classLabels = modules,
    identifiers = background,
    refCollection = GOcoll,
    useBackground = "given",
    nBestDataSets = length(GOcoll$dataSets),
    # threshold = 0.1,
    # thresholdType = "Bonferroni",
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = "not_significant",
    maxReportedOverlapGenes = 500
  )
  
  enrichmentTable <- GOenrichment$enrichmentTable %>%
    filter(nCommonGenes > 10, pValue <= 0.1)
  
  return(enrichmentTable)
}

GOcollection <- buildGOcollection(organism = "mouse")
# GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")
background <- fread(file = "data/kimono_input/prior_expr_funcoup_mm.csv")
background <- unique(c(background$Gene_A, background$Gene_B))
background <- mapIds(org.Mm.eg.db, keys = background, 
                     keytype = "ENSEMBL", column="ENTREZID")
for (reg in regions){
  
  go_enr_unique <- go_enrichment_all(list_reg[[reg]], GOcollection, TRUE, background)
  fwrite(go_enr_unique, file = paste0("tables/coExpression_kimono/03_AnalysisFuncoup/",
                                      "/09_", reg, "_GOterms_unique_", mode, ".csv"))
  go_enr_all <- go_enrichment_all(list_reg[[reg]], GOcollection, FALSE, background)
  fwrite(go_enr_all, file = paste0("tables/coExpression_kimono/03_AnalysisFuncoup/",
                                      "/09_", reg, "_GOterms_all_", mode, ".csv"))
  
  list_reg[[reg]]$GOterms_unique <- sapply(list_reg[[reg]]$entrez,
                                           function(x) paste(go_enr_unique$dataSetName[which(str_detect(go_enr_unique$overlapGenes, x))], collapse="|"))
  list_reg[[reg]]$GOterms_all <- sapply(list_reg[[reg]]$entrez,
                                        function(x) paste(go_enr_all$dataSetName[which(str_detect(go_enr_all$overlapGenes, x))], collapse="|"))
}


# 4. Print df of each brain region to file -------------------

for (reg in regions){
  # list_reg[[reg]]$ensembl_id <- NULL
  fwrite(list_reg[[reg]], file = paste0("tables/coExpression_kimono/03_AnalysisFuncoup/",
                                        "/09_", reg, "_hubGenes_unique_GOterms_", mode, ".csv"))
}

