##################################################
## Project: DexStim Mouse Brain
## Date: 21.09.2020
## Author: Nathalie
##################################################
# GO plots single regions

setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

library(org.Mm.eg.db)
library(dplyr)
library(stringr)
library(anRichment)
library(anRichmentMethods)

regions <- c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1")

folder_plots <- paste0("figures")
folder_tables <- paste0("tables")


# 1. read DE tables from all regions ----------

list_reg <- list()
for (reg in regions){
  res <- read.table(file=paste0(folder_tables, "/02_", reg, "_deseq2_Dex_1_vs_0_lfcShrink.txt"),sep="\t")
  res <- res[res$padj <= 0.1,]
  res$ensembl_id <- rownames(res)
  # res$padj[which(is.na(res$padj))] <- 1
  res$gene_symbol <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                            keytype = "ENSEMBL", column="SYMBOL")
  res$entrez <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                       keytype = "ENSEMBL", column="ENTREZID")
  list_reg[[reg]] <- res
}


# 2. check uniqueness of DE genes ---------------

for (reg in regions){
  index_reg <- which(regions == reg)
  # df <- bind_rows(list_reg[-index_reg], .id="region") %>%
  #   filter(DE0.1)
  df <- bind_rows(list_reg[-index_reg], .id="region")
  list_reg[[reg]]$regions_DE <- sapply(list_reg[[reg]]$ensembl_id, 
                                       function(x) paste(df[df$ensembl_id == x,]$region, collapse = " "))
  list_reg[[reg]]$unique_DE <- sapply(list_reg[[reg]]$regions_DE, 
                                      function(x) x == "")
}


# 3. GO enrichment for the genes of each region ------------------

go_enrichment_all <- function(df_reg, GOcoll, unique){
  if (unique){
    genes <- df_reg$entrez[df_reg$unique_DE]
  } else {
    genes <- df_reg$entrez
  }
  background <- read.table(file = paste0(folder_tables, "/06_background_entrezID.txt"),
                           header = FALSE)
  modules <- rep("not_significant", nrow(background))
  modules[which(background$V1 %in% genes)] <- "significant"
  
  # enrichment
  GOenrichment <- enrichmentAnalysis(
    classLabels = modules,
    identifiers = background$V1,
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
  
  enrichmentTable <- GOenrichment$enrichmentTable
  
  return(enrichmentTable)
}

GOcollection <- buildGOcollection(organism = "mouse")
list_GO <- list()
# GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")
for (reg in regions){
  
  go_enr_unique <- go_enrichment_all(list_reg[[reg]], GOcollection, TRUE)
  # go_enr_all <- go_enrichment_all(list_reg[[reg]], GOcollection, FALSE)
  list_GO[[reg]] <- go_enr_unique
  
}


# 4. Plot GO terms
df_all <- bind_rows(list_GO, .id="region")
for (reg in regions){
  
  df_reg <- list_GO[[reg]] %>%
    filter(nCommonGenes >= 10, pValue <= 0.1) %>%
    group_by(inGroups) %>% slice_min(order_by = pValue, n = 10)

  df <- df_all[df_all$dataSetName %in% df_reg$dataSetName,]
  df$dataSetName <- sapply(df$dataSetName, function(x) str_trunc(x, 45, "right"))
  df$dataSetName <- factor(df$dataSetName, levels = rev(reorder(df$dataSetName[df$region==reg], df$pValue[df$region==reg])))
  df$region <- factor(df$region, levels = c("AMY", "CER", "PFC", "PVN", "vDG", "dDG", "vCA1", "dCA1"))
  df$inGroups <- factor(df$inGroups)
  levels(df$inGroups) <- c("Biological Process", "Cellular Components", "Molecular Function")
  
  # Plot results (plotted pvalues are not adjusted for multiple testing)
  df %>%
    arrange(desc(dataSetName)) %>%
    ggplot(aes(x=dataSetName, y = -log10(pValue), fill = region)) +
    geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity") +
    geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
    coord_flip() +
    xlab("GOterm") +
    ggtitle(paste0("GO terms enriched for DE genes only in ",reg, " (Top 10 each)")) +
    facet_wrap(~inGroups, scales="free") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(filename = paste0(folder_plots, "/10_", reg, "_GOterms_unique.png"), width = 13, height = 7)
}


