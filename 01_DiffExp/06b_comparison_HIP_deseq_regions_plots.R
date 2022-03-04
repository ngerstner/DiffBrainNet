##################################################
## Project: DexStim Mouse Brain
## Date: 19.04.2021
## Author: Nathalie
##################################################
# Compare deseq SV DE genes between HIP regions
# Plots for manuscript

library(RColorBrewer)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
# library(gridExtra)
library(pheatmap)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"

regions <-
  c("vDG", "dDG", "vCA1", "dCA1")



# 1. Read DE tables from HIP regions ----------

list_reg <- list()
list_reg_sig <- list()
list_genes_top10 <- list()

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
  list_reg[[reg]] <- res
  res_sig <- res[res$padj <= 0.1, ]
  # res_sig <- res[res$log2FoldChange >= 1]
  list_reg_sig[[reg]] <- res_sig
  list_genes_top10[[reg]] <- res_sig$Ensembl_ID[1:25]
}


# 2. check uniqueness of DE genes ---------------

for (reg in regions){
  index_reg <- which(regions == reg)
  # df <- bind_rows(list_reg[-index_reg], .id="region") %>%
  #   filter(DE0.1)
  df <- bind_rows(list_reg_sig[-index_reg], .id="region")
  list_reg_sig[[reg]]$regions_DE <- sapply(list_reg_sig[[reg]]$Ensembl_ID, 
                                       function(x) paste(df[df$Ensembl_ID == x,]$region, collapse = " "))
  list_reg_sig[[reg]]$unique_DE <- sapply(list_reg_sig[[reg]]$regions_DE, 
                                      function(x) x == "")
}


# 3. plot number of unique/not unique genes within HIP -------------------------

# make dataframe to plot barplot
df_uni <- data.frame(region = character(),
                     uniqueness = character(),
                     genes = numeric())
for (reg in regions){
  u <- table(list_reg_sig[[reg]]$unique_DE)
  df_uni <- rbind(df_uni, list("region" = reg, "uniqueness" = "unique", 
                               "genes" = u[2]))
  df_uni <- rbind(df_uni, list("region" = reg, "uniqueness" = "not_unique", 
                               "genes" = u[1]))
}
df_uni$region <- as.factor(df_uni$region)
df_uni$uniqueness <- as.factor(df_uni$uniqueness)

# sum DE genes per region for percentage label
df_uni <- df_uni %>%
  group_by(region) %>%
  mutate(sum = sum(genes))

# stacked barplot with unique DE genes in HIP
ggplot(df_uni, aes(x = region, 
                   y = genes,
                   fill = uniqueness)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    name = "",
    labels = c("DE in multiple HIP regions", "DE unique in HIP"),
    values = c("tan1", "royalblue")
  ) +
  xlab("brain region") +
  ylab("# diff. exp. genes") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_blank(),
    legend.position = "top"
  ) +
  geom_text(aes(label = paste0(round((genes / sum) * 100, digits = 1
  ), "%")),
  position = position_stack(vjust = 0.5),
  size = 4)
ggsave(filename = paste0(basepath, "figures/06b_comparison_HIP_deseq_barplot_vertical_percentage.png"),
       width = 8,
       height = 6)



# 4. Heatmap of unique genes with highest p-value ---------------------------

# Function to get gene IDs of unique DE with lowest p-values
top_unique <- function(df){
  
  df <- df[df$unique_DE,]
  genes <- df$Ensembl_ID[1:10]
  
  return(genes)
}

# Concatenate DE tables -----------------
# TODO: try also p-value

data <- bind_rows(list_reg, .id = "region")
genes <- unlist(lapply(list_reg_sig, top_unique))
data <- data[data$Ensembl_ID %in% genes,] %>%
  pivot_wider(id_cols = c(Gene_Symbol),
              names_from = region,
              #values_from = log2FoldChange)
              values_from = padj)

data_mat <- as.matrix(data[,2:5])
rownames(data_mat) <- data$Gene_Symbol


# Heatmap 

# Complex heatmap
library(ComplexHeatmap)
library(circlize)
Heatmap(data_mat, 
        name = "log2FoldChange", #title of legend
        column_title = "Hippocampal region", row_title = "Genes",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        cluster_columns = FALSE,
        #col = colorRamp2(c(-4, 0, 4), c("blue", "#EEEEEE", "red"))
        col = colorRamp2(c(0,1), c("red", "#EEEEEE"))
) 

pheatmap(data_mat, 
         cutree_rows = 2,
         cluster_cols = FALSE)




# 3. GO enrichment for the genes of each region ------------------

go_enrichment_all <- function(df_reg, unique){
  if (unique){
    genes <- df_reg$Ensembl_ID[df_reg$unique_DE]
  } else {
    genes <- df_reg$Ensembl_ID
  }
  genes <- mapIds(org.Mm.eg.db, keys = genes, keytype = "ENSEMBL",
                  column = "ENTREZID")
  background <- read.table(file = paste0(basepath, "tables/06_background_entrezID.txt"),
                           header = FALSE)$V1
  
  # enrichment
  ego <- enrichGO(gene          = as.character(genes),
                  universe      = as.character(background),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  #pvalueCutoff  = 0.01,
                  #qvalueCutoff  = 0.05,
                  minGSSize = 10,    # min number of genes associated with GO term
                  maxGSSize = 10000, # max number of genes associated with GO term
                  readable      = TRUE)@result
  # ego_simple <- clusterProfiler::simplify(
  #   ego,
  #   cutoff = 0.7,
  #   by = "p.adjust",
  #   select_fun = min,
  #   measure = "Wang",
  #   semData = NULL
  # )@result
  
  return(ego)
}

list_GO <- list()
# GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")
for (reg in regions){
  
  go_enr_unique <- go_enrichment_all(list_reg_sig[[reg]], TRUE)
  # go_enr_all <- go_enrichment_all(list_reg[[reg]], GOcollection, FALSE)
  list_GO[[reg]] <- go_enr_unique
  
}


# 4. Plot GO terms
df_all <- bind_rows(list_GO, .id="region")
for (reg in regions){
  
  df_reg <- list_GO[[reg]] %>%
    # group_by(inGroups) %>% 
    # slice_min(order_by = p.adjust, n = 10)
    slice_head(n = 10)
  
  df <- df_all[df_all$ID %in% df_reg$ID,]
  # df$dataSetName <- sapply(df$dataSetName, function(x) str_trunc(x, 45, "right"))
  df$Description <- factor(df$Description, 
                           levels = rev(reorder(df$Description[df$region==reg], df$p.adjust[df$region==reg])))
  df$region <- factor(df$region, levels = c("vDG", "dDG", "vCA1", "dCA1"))
  #df$inGroups <- factor(df$inGroups)
  #levels(df$inGroups) <- c("Biological Process", "Cellular Components", "Molecular Function")
  
  # Plot results (plotted pvalues are not adjusted for multiple testing)
  df %>%
    arrange(desc(Description)) %>%
    ggplot(aes(x=Description, y = -log10(p.adjust), fill = region)) +
    geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity") +
    geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
    coord_flip() +
    xlab("GOterm") +
    ggtitle(paste0("GO terms enriched for DE genes only in ",reg, " (Top 10 each)")) +
    # facet_wrap(~inGroups, scales="free") +
    theme(axis.text.y = element_text(size = 10))
  ggsave(filename = paste0(basepath, "figures", "/06b_HIP_", reg, "_GOterms_unique.png"), width = 13, height = 7)
}


