##################################################
## Project: DexStim Mouse Brain
## Date: 05.10.2021
## Author: Nathalie
##################################################
# Compare GO enrichment of unique DE and hub genes

library(data.table)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
regions <- c("AMY", "CER", "dCA1", "dDG", "PFC", "PVN", "vCA1", "vDG")
beta_cutoff <- 0.01

folder_plots <- paste0("figures")
folder_tables <- paste0("tables")


# 1. Read DE and hub genes from all regions and background ----------

# DE genes
list_reg_de <- list()
for (reg in regions){
  res <- read.table(file=paste0(basepath, folder_tables, "/02_", reg, "_deseq2_Dex_1_vs_0_lfcShrink.txt"),
                    sep="\t", header = TRUE)
  res <- res[res$padj <= 0.1,]
  res$gene_symbol <- mapIds(org.Mm.eg.db, keys = res$Ensembl_ID, 
                            keytype = "ENSEMBL", column="SYMBOL")
  res$entrez <- mapIds(org.Mm.eg.db, keys = res$Ensembl_ID, 
                       keytype = "ENSEMBL", column="ENTREZID")
  list_reg_de[[reg]] <- res
}

# hub genes
list_reg_hub <- list()
for (reg in regions){
  res <- fread(paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/04_",
                      reg,"_funcoup_differential_nodeBetweennessNorm_betacutoff",beta_cutoff,".csv"))
  res <- res[res$nodebetweenness_norm>=1.0 & ! is.na(res$nodebetweenness_norm)]
  res$entrez <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                       keytype = "ENSEMBL", column="ENTREZID")
  list_reg_hub[[reg]] <- res
}

# background are all genes in out dataset
background <- read.table(file = paste0(basepath, "tables/06_background_entrezID.txt"),
                         header = FALSE)[,1]


# 2. check uniqueness of DE and hub genes ---------------

# DE genes
for (reg in regions){
  index_reg <- which(regions == reg)
  df <- bind_rows(list_reg_de[-index_reg], .id="region")
  list_reg_de[[reg]]$regions_DE <- sapply(list_reg_de[[reg]]$Ensembl_ID, 
                                       function(x) paste(df[df$Ensembl_ID == x,]$region, collapse = " "))
  list_reg_de[[reg]]$unique_DE <- sapply(list_reg_de[[reg]]$regions_DE, 
                                      function(x) x == "")
}

# hub genes
for (reg in regions){
  index_reg <- which(regions == reg)
  df <- bind_rows(list_reg_hub[-index_reg], .id="region")
  list_reg_hub[[reg]]$regions_hub <- sapply(list_reg_hub[[reg]]$ensembl_id,
                                          function(x) paste(df[df$ensembl_id == x,]$region, collapse = " "))
  list_reg_hub[[reg]]$unique_hub <- sapply(list_reg_hub[[reg]]$regions_hub, 
                                         function(x) x == "")
}


# 2. Plot top GO enrichment of DE and hub genes per brain regions

#minCount <- 0
#minCount <- 5
# minCount <- 10

# --> changed to minCount with regard to number of genes in geneset

for (reg in regions){
  
  # Unique DE genes
  genes <- list_reg_de[[reg]]$entrez[list_reg_de[[reg]]$unique_DE]
  
  # GO enrichment for unique DE genes
  ego_de <- enrichGO(gene          = as.character(genes),
                  universe      = as.character(background),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  # pvalueCutoff  = 1,
                  # qvalueCutoff  = 1,
                  minGSSize = 5,    # min number of genes associated with GO term
                  maxGSSize = 10000, # max number of genes associated with GO term
                  readable      = TRUE)@result
  # min number of genes overlapping
  min_count <- ceiling(length(genes)*0.15)
  ego_de_filt <- ego_de[ego_de$Count >= min_count,]
  
  
  # Unique hub genes
  genes <- list_reg_hub[[reg]]$entrez[list_reg_hub[[reg]]$unique_hub]
  
  # GO enrichment for unique hub genes
  ego_hub <- enrichGO(gene          = as.character(genes),
                     universe      = as.character(background),
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     # pvalueCutoff  = 1,
                     # qvalueCutoff  = 1,
                     minGSSize = 5,    # min number of genes associated with GO term
                     maxGSSize = 10000, # max number of genes associated with GO term
                     readable      = TRUE)@result
  # min number of genes overlapping
  min_count <- ceiling(length(genes)*0.15)
  ego_hub_filt <- ego_hub[ego_hub$Count >= min_count,]
  
  
  # Plot top 10 genes of DE and hub with comparison to each other
  
  # Top 10 terms for DE genes
  ego_de_top10 <- head(ego_de_filt, n = 10) %>%
    dplyr::mutate("main" = TRUE)
  ego_de_top10_hub <- ego_hub[ego_hub$ID %in% ego_de_top10$ID,] %>%
    dplyr::mutate("main" = FALSE)
  
  ego_de_plot <- bind_rows("de" = ego_de_top10, "hub" = ego_de_top10_hub,
                           .id = "groups")
  
  # Top 10 terms for hub genes
  ego_hub_top10 <- head(ego_hub_filt, n = 10) %>%
    dplyr::mutate("main" = TRUE)
  ego_hub_top10_de <- ego_de[ego_de$ID %in% ego_hub_top10$ID,] %>%
    dplyr::mutate("main" = FALSE)
  
  ego_hub_plot <- bind_rows("hub" = ego_hub_top10, "de" = ego_hub_top10_de,
                            .id = "groups")
  
  # Combine them into one dataframe
  ego_plot <- bind_rows("de" = ego_de_plot, "hub" = ego_hub_plot,
                        .id = "subplot")
  ego_plot$Description <- factor(ego_plot$Description, 
                                 levels = unique(ego_plot$Description))
  
  # labeller function for facet titles
  facet_names <- c(
    'de'="Top 10 GO terms for DE genes",
    'hub'="Top 10 GO terms for hub genes"
  )
  
  # Plot them all together
  print(ggplot(ego_plot, aes(x=Description, y = -log10(p.adjust), 
                             fill = groups, colour = main)) +
          geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity") +
          geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
          scale_colour_manual(values = c("grey", "black"), guide = FALSE) +
          scale_fill_manual(name = "Geneset",
                            values = c("orange", "darkred"),
                            labels = c("DE genes", "hub genes")) +
          coord_flip() +
          scale_x_discrete(limits = rev) +
          xlab("GOterm") +
          ggtitle(paste0("GO terms enriched for DE and hub genes only in ",reg)) +
          facet_wrap(~subplot, scales="free", labeller = as_labeller(facet_names)) +
          theme(axis.text.y = element_text(size = 10)))
  # ggsave(filename = paste0(basepath, folder_plots, "/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
  #                          "10_",reg,"_GOterms_DEandHubGenes_min", minCount,".png"),
  #        width = 13, height = 7)
  ggsave(filename = paste0(basepath, folder_plots, "/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                           "10_",reg,"_GOterms_DEandHubGenes.pdf"),
         width = 11, height = 7)
  
  
  
}


