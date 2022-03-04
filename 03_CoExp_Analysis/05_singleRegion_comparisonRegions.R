##################################################
## Project: DexStim Mouse Brain
## Date: 15.12.2020
## Author: Nathalie
##################################################
# Use beta cutoff and analyze top genes (nodebetweenness)

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(eulerr)
library(UpSetR)
library(org.Mm.eg.db)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
regions <- c("AMY", "CER", "dCA1", "dDG", "PFC", "PVN", "vCA1", "vDG")
beta_cutoff <- 0.01


# 0. functions -------------------------------
write_genelist <- function(genelist, filename){
  # write list with ENSEMBL IDs
  write.table(genelist, file = paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                                      "/05_",filename,"_ensemblID.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  # write list with ENTREZ IDs
  entrez <- mapIds(org.Mm.eg.db, keys = genelist, keytype = "ENSEMBL", column="ENTREZID")
  write.table(entrez, file = paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                                    "/05_",filename,"_entrezID.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  # write list with GENE SYMBOLS
  symbol <- mapIds(org.Mm.eg.db, keys = genelist, keytype = "ENSEMBL", column="SYMBOL")
  write.table(symbol, file = paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
                                    "/05_",filename,"_geneSymbol.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# 1. Read data --------------------------------
nodedegrees_list <- list()
nodedegrees_0.5 <- list()
for (reg in regions){
  nodedegrees_list[[reg]] <- fread(paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                        "04_",reg,"_funcoup_differential_nodedegreesNorm_betacutoff",beta_cutoff,".csv"))
  nodedegrees_0.5[[reg]] <- nodedegrees_list[[reg]]$ensembl_id[nodedegrees_list[[reg]]$nodedegree_norm>=0.5 &
                                                                 ! is.na(nodedegrees_list[[reg]]$nodedegree_norm)]
  
}

nodebetweenness_list <- list()
nodebetweenness_1 <- list()
for (reg in regions){
  nodebetweenness_list[[reg]] <- fread(paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                        "04_",reg,"_funcoup_differential_nodebetweennessNorm_betacutoff",beta_cutoff,".csv"))
  nodebetweenness_1[[reg]] <- nodebetweenness_list[[reg]]$ensembl_id[nodebetweenness_list[[reg]]$nodebetweenness_norm>=1 &
                                                                       ! is.na(nodebetweenness_list[[reg]]$nodebetweenness_norm)]
}


# 2. Compare top genes between regions using Upset Plot

# 2.1 Nodedegree
png(filename = paste0(basepath, "/figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                      "upsetPlot_normNodedegree0.5.png"),
    height = 700, width = 1000)
print(upset(fromList(nodedegrees_0.5), nsets = 8, nintersects = 50, order.by = "freq",
            text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8)))
dev.off()

# 2.2 Nodebetweenness
png(filename = paste0(basepath, "/figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                      "upsetPlot_normNodebetweenness1.0.png"),
    height = 700, width = 1000)
print(upset(fromList(nodebetweenness_1), nsets = 8, nintersects = 50, order.by = "freq",
            text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8)))
dev.off()



# 3. Plot correlation between 2 different brain regions

reg_comb <- combn(regions, 2)
# reg1 <- "PFC"
# reg2 <- "dDG"

for (i in 1:ncol(reg_comb)){
  
  reg1 <- reg_comb[1,i]
  reg2 <- reg_comb[2,i]
  
  # 3.1 Nodedegree
  degree_reg <- inner_join(nodedegrees_list[[reg1]], nodedegrees_list[[reg2]], by = "ensembl_id",
                           suffix = c(".reg1", ".reg2"))
  ggplot(degree_reg, aes(x=nodedegree_norm.reg1, y=nodedegree_norm.reg2)) +
    # geom_point(size=1,alpha = 0.1)
    geom_hex() +
    # geom_bin2d() +
    xlab(paste("norm. nodedegree", reg1)) +
    ylab(paste("norm. nodedegree", reg2)) +
    ggtitle(paste("Nodedegree in", reg1, "and", reg2, "differential network"))
  ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                           "comparison_",reg1,"-",reg2,"_correlationNodedegreeNorm.png"),
         width = 9, height = 8)
  
  # 3.2 Nodebetweenness
  between_reg <- inner_join(nodebetweenness_list[[reg1]], nodebetweenness_list[[reg2]], by = "ensembl_id",
                            suffix = c(".reg1", ".reg2"))
  ggplot(between_reg, aes(x=nodebetweenness_norm.reg1, y=nodebetweenness_norm.reg2)) +
    # geom_point(size=1,alpha = 0.1)
    geom_hex() +
    # geom_bin2d() +
    xlab(paste("norm. nodebetweenness", reg1)) +
    ylab(paste("norm. nodebetweenness", reg2)) +
    ggtitle(paste("Nodebetweenness in", reg1, "and", reg2, "differential network"))
  ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                           "comparison_",reg1,"-",reg2,"_correlationNodebetweennessNorm.png"),
         width = 9, height = 8)
  
}


# 4. Gene lists -------------------------------------

# 4.1 Overlap all regions
# nodedegree
overlap_degree <- Reduce(intersect, nodedegrees_0.5)
write_genelist(overlap_degree, "topgenesNodedegree0.5_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1")
# nodebetweenness
overlap_between <- Reduce(intersect, nodebetweenness_1)
write_genelist(overlap_degree, "topgenesNodebetweenness1_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1")

# union of all nodebetweenness hub genes
union_between <- Reduce(union, nodebetweenness_1)

# comparison with union of all DE genes
de_genes <- fread(paste0(basepath, "tables/06_union_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.txt"),
                  header = FALSE)

# overlap of hub genes and de genes
common_hub_de <- intersect(union_between, de_genes$V1)
