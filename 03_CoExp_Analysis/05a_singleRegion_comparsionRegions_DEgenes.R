##################################################
## Project: DexStim Mouse Brain
## Date: 15.12.2020
## Author: Nathalie
##################################################
# Use beta cutoff and analyze top genes (nodebetweenness)
# UpSet Plot of DE genes with nodebetweenness >/< 1

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


# 1. Read data --------------------------------
# nodedegrees_list <- list()
# nodedegrees_0.5 <- list()
# for (reg in regions){
#   nodedegrees_list[[reg]] <- fread(paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
#                                           "04_",reg,"_funcoup_differential_nodedegreesNorm_betacutoff",beta_cutoff,".csv"))
#   nodedegrees_0.5[[reg]] <- nodedegrees_list[[reg]]$ensembl_id[nodedegrees_list[[reg]]$nodedegree_norm>=0.2 &
#                                                                  ! is.na(nodedegrees_list[[reg]]$nodedegree_norm)]
#   
# }

de_nodebetween <- list()
for (reg in regions){
  nodebetweenness <- fread(paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                              "04_",reg,"_funcoup_differential_nodebetweennessNorm_betacutoff",beta_cutoff,".csv"))
  nodebetweenness_1 <- nodebetweenness$ensembl_id[nodebetweenness$nodebetweenness_norm>= 1.0&
                                                                       ! is.na(nodebetweenness$nodebetweenness_norm)]
  
  de_genes <- fread(paste0(basepath, "tables/02_",reg,"_deseq2_Dex_1_vs_0_lfcShrink.txt")) %>%
    filter(padj <= 0.1)
  
  de_nodebetween[[reg]] <- intersect(nodebetweenness_1, de_genes$Ensembl_ID)
}


# 2 Upset Plot
png(filename = paste0(basepath, "/figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                      "upsetPlot_DE_normNodebetweennessAbove1.0.png"),
    height = 700, width = 1000)
print(upset(fromList(de_nodebetween), nsets = 8, nintersects = 50, order.by = "freq",
            text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8)))
dev.off()




# 1b. Read data --------------------------------

de_nodebetween <- list()
for (reg in regions){
  nodebetweenness <- fread(paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                  "04_",reg,"_funcoup_differential_nodebetweennessNorm_betacutoff",beta_cutoff,".csv"))
  nodebetweenness_1 <- nodebetweenness$ensembl_id[nodebetweenness$nodebetweenness_norm>= 1.0&
                                                    ! is.na(nodebetweenness$nodebetweenness_norm)]
  
  de_genes <- fread(paste0(basepath, "tables/02_",reg,"_deseq2_Dex_1_vs_0_lfcShrink.txt")) %>%
    filter(padj <= 0.1)
  
  de_nodebetween[[reg]] <- nodebetweenness_1
  de_nodebetween[[paste0(reg,"_de")]] <- de_genes$Ensembl_ID
}

df <- fromList(de_nodebetween)
# x <- which(df$PFC == 1 & rowSums(df[,seq(1,15,by=2)]) == 1)
# df$de <- rowSums(df[,seq(2,16,by=2)])
# y <- which(df$de >= 1)

# function to count genes that are also DE gene in at least one 
# of the intersect regions
de_region <- function(x){
  index_de <- which(x[seq(1,15,by=2)] == 1)*2
  s <- sum(x[index_de])
  return(s)
}
df$de_region <- apply(df, 1, de_region)

# 2b Upset Plot
pdf(file = paste0(basepath, "/figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                      "upsetPlot_DEcolor_normNodebetweennessAbove1.0.pdf"),
    height = 10, width = 14)
print(upset(df, nsets = 16, nintersects = 50, order.by = "freq", 
            sets = c("AMY", "CER", "dCA1",
                     "dDG", "PFC", "PVN",
                     "vCA1", "vDG"),
            text.scale = c(1.8, 1.9, 1.8, 1.9, 1.9, 1.9),
            sets.x.label = "#Hub genes in brain region",
            mainbar.y.label = "#Hub genes in intersection",
            queries = list(
              list(
                query = elements,
                params = list("de_region",1),
                color = "#FFA500", 
                active = T,
                query.name = "DE gene in at least one of the intersect regions"
              )
            ),
            query.legend = "bottom"))
dev.off()



# 2. Read nodedegree data --------------------------------

de_nodedegree <- list()
for (reg in regions){
  nodedegree <- fread(paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                  "04_",reg,"_funcoup_differential_nodedegreesNorm_betacutoff",beta_cutoff,".csv"))
  nodedegree_1 <- nodedegree$ensembl_id[nodedegree$nodedegree_norm>= 0.5&
                                                    ! is.na(nodedegree$nodedegree_norm)]
  
  de_genes <- fread(paste0(basepath, "tables/02_",reg,"_deseq2_Dex_1_vs_0_lfcShrink.txt")) %>%
    filter(padj <= 0.1)
  
  de_nodedegree[[reg]] <- nodedegree_1
  de_nodedegree[[paste0(reg,"_de")]] <- de_genes$Ensembl_ID
}

df <- fromList(de_nodedegree)
# x <- which(df$PFC == 1 & rowSums(df[,seq(1,15,by=2)]) == 1)
# df$de <- rowSums(df[,seq(2,16,by=2)])
# y <- which(df$de >= 1)

# function to count genes that are also DE gene in at least one 
# of the intersect regions
de_region <- function(x){
  index_de <- which(x[seq(1,15,by=2)] == 1)*2
  s <- sum(x[index_de])
  return(s)
}
df$de_region <- apply(df, 1, de_region)

# 2b Upset Plot
png(filename = paste0(basepath, "/figures/02_CoExp_Kimono/03_AnalysisFuncoup/comparisonRegions/",
                      "upsetPlot_DEcolor_normNodedegreeAbove0.5.png"),
    height = 700, width = 1000)
print(upset(df, nsets = 8, nintersects = 50, order.by = "freq", 
            sets = c("AMY", "CER", "dCA1",
                     "dDG", "PFC", "PVN",
                     "vCA1", "vDG"),
            text.scale = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.8),
            queries = list(
              list(
                query = elements,
                params = list("de_region",1),
                color = "#Df5286", 
                active = T,
                query.name = "DE gene in at least one of the intersect regions"
              )
            ),
            query.legend = "bottom"))
dev.off()

