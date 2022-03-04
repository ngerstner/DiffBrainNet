##################################################
## Project: DexStim Mouse Brain
## Date: 08.01.2020
## Author: Nathalie
##################################################
# Comparison of important and region specific genes
# according to DE and network analysis (nodebetweenness)

library(org.Mm.eg.db)
library(data.table)
library(ggplot2)
library(dplyr)
library(eulerr)
library(gridExtra)
library(grid)
library(igraph)
library(RCy3)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
regions <- c("AMY", "CER", "dCA1", "dDG", "PFC", "PVN", "vCA1", "vDG")
beta_cutoff <- 0.01


### ANALYSIS -----------------------------

# 1. Read data from all regions ----------
# 1a. DE tables

list_de <- list()
for (reg in regions){
  res <- read.table(file=paste0(basepath, "tables/02_", reg, 
                                "_deseq2_Dex_1_vs_0_lfcShrink.txt"),sep="\t")
  res <- res[res$padj <= 0.1,]
  res$ensembl_id <- rownames(res)
  # res$padj[which(is.na(res$padj))] <- 1
  res$gene_symbol <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
                            keytype = "ENSEMBL", column="SYMBOL")
  list_de[[reg]] <- res
}

# # 1b. Network tables
# list_net <- list()
# for (reg in regions){
#   res <- fread(paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
#                        "04_", reg, "_funcoup_differential_nodebetweennessNorm_betacutoff", 
#                        beta_cutoff, ".csv"))
#   res <- res[res$nodebetweenness_norm>=1.0 & ! is.na(res$nodebetweenness_norm)]
#   list_net[[reg]] <- res
# }
# 
# 
# 
# # 2. Venn/Euler Plot per region ---------------
# list_euler <- list()
# 
# for (reg in regions){
# 
#   list1 <- list(list_de[[reg]]$ensembl_id,
#                 list_net[[reg]]$ensembl_id)
#   names(list1) <- c("Diff. exp. genes",
#                     "Diff. co-exp. genes")
#   list_euler[[reg]] <- plot(euler(list1, shape = "ellipse"),
#              labels = list(cex = 1.0), quantities = list(cex = 1.0),
#              main = paste0(reg))
# 
# }
# 
# grid.arrange(grobs = list_euler, ncol = 4,
#              top = "Comparison of DE genes and top network genes")





### ANALYZE NEIGHBOURS OF DE GENES IN NETWORK #######################

# 1b. Network tables
list_diff <- list()
for (reg in regions){
  res <- fread(paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                      "04_singleRegion_", reg, "_filtered_diffNetwork.csv"))
  # res <- res[res$nodebetweenness_norm>=1.0 & ! is.na(res$nodebetweenness_norm)]
  # res$entrez <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
  #                      keytype = "ENSEMBL", column="ENTREZID")
  list_diff[[reg]] <- res
}

list_base <- list()
for (reg in regions){
  res <- fread(paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                      "04_singleRegion_", reg, "_filtered_baselineNetwork.csv"))
  # res <- res[res$nodebetweenness_norm>=1.0 & ! is.na(res$nodebetweenness_norm)]
  # res$entrez <- mapIds(org.Mm.eg.db, keys = res$ensembl_id, 
  #                      keytype = "ENSEMBL", column="ENTREZID")
  list_base[[reg]] <- res
}

# 1c. Network nodebetweeness tables
list_net_base <- list()
list_net_diff <- list()
for (reg in regions){
  res <- fread(paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                      "04_", reg, "_funcoup_differential_nodebetweennessNorm_betacutoff", beta_cutoff, ".csv"))
  list_net_diff[[reg]] <- res
  res <- fread(paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                      "03_", reg, "_funcoup_baseline_nodebetweennessNorm_betacutoff", beta_cutoff, ".csv"))
  list_net_base[[reg]] <- res
}



# 2. Subset networks to DE genes 
for (reg in regions){
  # norm_nodebetweenness is NA whenever prior nodebetweenness < 10000
  # (our definition)
  de_genes <- data.frame("ensembl_id" = list_de[[reg]]$ensembl_id) %>%
    left_join(list_net_diff[[reg]], by = "ensembl_id", ) %>%
    left_join(list_net_base[[reg]], by = c("ensembl_id", "gene_symbol"), 
              suffix = c(".diff", ".base"))
  de_genes <- as_tibble(de_genes)
  
  # identify baseline neighbours of de_genes using igraph
  actors <- data.frame(name=unique(c(list_base[[reg]]$target, 
                                     list_base[[reg]]$predictor,
                                     de_genes$ensembl_id)))
  relations <- data.frame(from=list_base[[reg]]$target,
                          to=list_base[[reg]]$predictor)
  ig <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  ig <- simplify(ig)
  # add baseline neighbors in column
  de_genes$neighbors.base <- sapply(de_genes$ensembl_id,
                               function(x) names(unlist(neighbors(ig, x))) )
  de_genes$nr_neigh_de.base <- sapply(de_genes$neighbors.base, 
                                 function(x) length(intersect(x, de_genes$ensembl_id)))
  de_genes$nr_neigh_notde.base <- sapply(de_genes$neighbors.base,
                                    function(x) length(setdiff(x, de_genes$ensembl_id)))
  de_genes$neighbors.base <- sapply(de_genes$neighbors.base, function(x) if(!length(x) == 0) 
    mapIds(org.Mm.eg.db, keys = x, keytype = "ENSEMBL", column = "SYMBOL"))
  de_genes$neighbors.base <- sapply(de_genes$neighbors.base, function(x) toString(x))
  
  
  # identify differential neighbours of de_genes using igraph
  actors <- data.frame(name=unique(c(list_diff[[reg]]$target, 
                                     list_diff[[reg]]$predictor,
                                     de_genes$ensembl_id)))
  relations <- data.frame(from=list_diff[[reg]]$target,
                          to=list_diff[[reg]]$predictor)
  ig <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  ig <- simplify(ig)
  # add differential neighbors in column
  de_genes$neighbors.diff <- sapply(de_genes$ensembl_id,
                                    function(x) names(unlist(neighbors(ig, x))) )
  de_genes$nr_neigh_de.diff <- sapply(de_genes$neighbors.diff, 
                                      function(x) length(intersect(x, de_genes$ensembl_id)))
  de_genes$nr_neigh_notde.diff <- sapply(de_genes$neighbors.diff,
                                         function(x) length(setdiff(x, de_genes$ensembl_id)))
  de_genes$neighbors.diff <- sapply(de_genes$neighbors.diff, function(x) if(!length(x) == 0) 
    mapIds(org.Mm.eg.db, keys = x, keytype = "ENSEMBL", column = "SYMBOL"))
  de_genes$neighbors.diff <- sapply(de_genes$neighbors.diff, function(x) toString(x))
  
  # move gene symbols to second column
  de_genes <- de_genes %>%
    select(ensembl_id, gene_symbol, everything())
  
  # write to file
  fwrite(de_genes, file = paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                                 "07_", reg, "_neighbors_DEgenes.csv"))
}