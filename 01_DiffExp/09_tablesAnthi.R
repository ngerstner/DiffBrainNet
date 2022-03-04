##################################################
## Project: DexStim Mouse Brain
## Date: 19.09.2020
## Author: Nathalie
##################################################
# Make tables for Anthi

setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

library(org.Mm.eg.db)
library(dplyr)
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
  # find regions where gene is also differentially expressed
  list_reg[[reg]]$regions_DE <- sapply(list_reg[[reg]]$ensembl_id, 
                                       function(x) paste(df[df$ensembl_id == x,]$region, collapse = " "))
  # boolean if gene is DE uniquely in this region
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
  
  enrichmentTable <- GOenrichment$enrichmentTable %>%
    filter(nCommonGenes > 10, pValue <= 0.1)
  
  return(enrichmentTable)
}

GOcollection <- buildGOcollection(organism = "mouse")
# GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")
for (reg in regions){
  
  go_enr_unique <- go_enrichment_all(list_reg[[reg]], GOcollection, TRUE)
  write.table(go_enr_unique, file = paste0(folder_tables, "/09_", reg, "_GOterms_unique.txt"),
            quote = FALSE, sep = "\t")
  go_enr_all <- go_enrichment_all(list_reg[[reg]], GOcollection, FALSE)
  write.table(go_enr_all, file = paste0(folder_tables, "/09_", reg, "_GOterms_all.txt"),
            quote = FALSE, sep = "\t")
  
  list_reg[[reg]]$GOterms_unique <- sapply(list_reg[[reg]]$entrez,
                                    function(x) paste(go_enr_unique$dataSetName[which(str_detect(go_enr_unique$overlapGenes, x))], collapse="|"))
  list_reg[[reg]]$GOterms_all <- sapply(list_reg[[reg]]$entrez,
                                           function(x) paste(go_enr_all$dataSetName[which(str_detect(go_enr_all$overlapGenes, x))], collapse="|"))
}


# 4. Print df of each brain region to file -------------------

for (reg in regions){
  list_reg[[reg]]$ensembl_id <- NULL
  write.csv(list_reg[[reg]], file = paste0(folder_tables, "/09_", reg, "_DEgenes_unique_GOterms.csv"),
              quote = FALSE)
}


# 5. Print logfoldchange in each region (examples for slides) -----------------------

# read all regions with all genes (no pval filtering)
list_reg <- list()
for (reg in regions){
  res <- read.table(file=paste0(folder_tables, "/02_", reg, "_deseq2_Dex_1_vs_0_lfcShrink.txt"),sep="\t")
  res$ensembl_id <- rownames(res)
  res$gene_symbol <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                            keytype = "ENSEMBL", column="SYMBOL")
  res$entrez <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                       keytype = "ENSEMBL", column="ENTREZID")
  list_reg[[reg]] <- res
}

# combine data of all regions (append rows)
df <- bind_rows(list_reg, .id="region")
head(df)
# df <- df %>%
#   dplyr::select(region, log2FoldChange, padj, ensembl_id, gene_symbol, entrez)

# all regions
genes_all <- read.table(file=paste0(folder_tables,"/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_symbolID.txt"))
for (i in 1:10){
  print(df %>%
  filter(gene_symbol == genes_all$V1[i]))
}
fkbp5 <- df %>%
  filter(gene_symbol == "Fkbp5")
ggplot(fkbp5, aes(x = region, y = log2FoldChange, fill = padj < 0.1)) +
  geom_bar(stat="identity") +
  xlab("brain region") +
  scale_fill_manual("FDR < 0.1", values = c("TRUE" = "yellowgreen", "FALSE" = "orange")) + 
  ggtitle("FKBP5: differentially expressed in all brain regions") +
  theme_bw() +
  theme(text = element_text(size=12))
ggsave(filename = paste0(folder_plots,"/09_FKBP5_foldchanges.png"), width = 6, height = 4)

# only CER
genes_cer <- read.table(file=paste0(folder_tables,"/06_unique_CER_symbolID.txt"))
for (i in 1:10){
  print(df %>%
          filter(gene_symbol == genes_cer$V1[i]))
}
tgfb3 <- df %>%
  filter(gene_symbol == "Tgfb3")
ggplot(tgfb3, aes(x = region, y = log2FoldChange, fill = padj < 0.1)) +
  geom_bar(stat="identity") +
  xlab("brain region") +
  scale_fill_manual("FDR < 0.1", values = c("TRUE" = "yellowgreen", "FALSE" = "orange")) + 
  ggtitle("TGFB3: differentially expressed only in the Cerebellum") +
  theme_bw() +
  theme(text = element_text(size=12))
ggsave(filename = paste0(folder_plots,"/09_TGFB3_foldchanges.png"), width = 6, height = 4)
