##################################################
## Project: DexStim Mouse Brain
## Date: 04.01.2020
## Author: Nathalie
##################################################
# Analyze top genes in baseline network
# Use beta cutoff and analyze top genes (nodebetweenness)

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(eulerr)
library(org.Mm.eg.db)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
region <- "PFC"
beta_cutoff <- 0.01
padj_cutoff <- 0.01
rsquared_cutoff <- 0.1

### FUNCTIONS -------------------------------------

# function to read all files from list
readFiles_concat <- function(file_list){
  
  # initialize empty data frame
  dataset <- data.frame()
  
  # read each file from list and append to data frame
  for (i in 1:length(file_list)){
    temp_data <- fread(file_list[i])
    dataset <- rbindlist(list(dataset, temp_data), use.names = T)
  }
  
  return(dataset)
}

# Z-score (z_ij) for the differential analysis between gene i and j
z_score <- function(beta_t, beta_c, se_t, se_c){
  z <- (beta_t - beta_c)/
    sqrt((se_t)^2 + (se_c)^2)
}

# Plot changes in ranks of nodes with highest betweenness
plotRanks <- function(a, b, labels = TRUE, labels.offset=0.1, arrow.len=0.1)
{
  old.par <- par(mar=c(1,1,1,1))
  
  # Find the length of the vectors
  len.1 <- length(a)
  len.2 <- length(b)
  
  # Plot two columns of equidistant points
  plot(rep(1, len.1), 1:len.1, pch=20, cex=0.8, 
       xlim=c(0, 3), ylim=c(0, max(len.1, len.2)),
       axes=F, xlab="", ylab="") # Remove axes and labels
  points(rep(2, len.2), 1:len.2, pch=20, cex=0.8)
  
  # Put labels next to each observation
  if (labels){
    text(rep(1-labels.offset, len.1), 1:len.1, a)
    text(rep(2+labels.offset, len.2), 1:len.2, b)
  }
  
  # Now we need to map where the elements of a are in b
  # We use the match function for this job
  a.to.b <- match(a, b)
  
  # Now we can draw arrows from the first column to the second
  arrows(rep(1.02, len.1), 1:len.1, rep(1.98, len.2), a.to.b, 
         length=arrow.len, angle=20)
  par(old.par)
}


### ANALYSIS ---------------------------------------

# 1. Prior and network
funcoup_prior <- fread(file = paste0(basepath, "data/kimono_input/prior_expr_funcoup_mm.csv"))
nodes_prior <- unique(c(funcoup_prior$Gene_A, funcoup_prior$Gene_B))
relations_prior <- data.frame(from = funcoup_prior$Gene_A,
                              to = funcoup_prior$Gene_B)
g_prior <- graph_from_data_frame(relations_prior, directed=FALSE, vertices=nodes_prior)
# Calculate nodedegrees of prior
nodedegree_prior <- sort(igraph::degree(g_prior), decreasing = TRUE)
# Calculate nodebetweenness of prior
# nodebetweenness_prior <- betweenness(g_prior, directed = FALSE)  # node betweenness: number of shortest paths going through a node
# saveRDS(nodebetweenness_prior, file = paste0(basepath, "data/workspaces/nodebetweenness_funcoup.rds"))
nodebetweenness_prior <- sort(readRDS(paste0(basepath, "data/workspaces/nodebetweenness_funcoup.rds")), decreasing = TRUE)
top100genes_prior <- names(nodebetweenness_prior[1:100])


# 2. Read kimono baseline expression networks
dex0_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_singleRegion\\_",region,"\\_dex0\\_funcoup\\_SE\\_.*\\.csv"),
                         full.names = TRUE)

data_dex0 <- readFiles_concat(dex0_files)


# 3. Remove interactions with very low r squared values & intercept & SVs
data <- data_dex0 %>% 
  filter(overall_rsq >= rsquared_cutoff) %>%
  filter(predictor != '(Intercept)')
data <- data[!startsWith(data$predictor, "SV"),]
# remove duplicated interactions (mistake made when separating nodes into chunks)
data <- data %>%
  distinct(target, predictor, .keep_all = TRUE)


# 4. Keep only interactions that have a beta value > cutoff
data <- data %>%
  mutate(betacut = (abs(beta_mean) > beta_cutoff))
data_cut <- data %>%
  filter(betacut)


# 5. Create baseline network corresponding to beta cutoff
head(data_cut[,c("target", "predictor", "beta_mean", "beta_stderr")], 20)
# Save filtered network
fwrite(data_cut, file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                "04_singleRegion_",region,"_filtered_baselineNetwork.csv"))

# Nodes in network
node_vec <- unique(c(data_cut$target, data_cut$predictor))

# Find modules
# relations <- data.frame(from=data_cut$target,
#                         to=data_cut$predictor,
#                         value=data_cut$z,
#                         performance=data_cut$p_adj)
relations <- data.frame(from=data_cut$target,
                        to=data_cut$predictor)
g_base <- graph_from_data_frame(relations, directed=FALSE, vertices=node_vec)
# does graph contain multiple edges with same start and endpoint or loop edges
is_simple(g_base)
g_base <- simplify(g_base) # check the edge attribute parameter

# Calculate nodedegree
nodedegree <- igraph::degree(g_base)
nodedegree <- sort(nodedegree, decreasing = TRUE)

# Calculate nodebetweenness
nodebetweenness <- betweenness(g_base, directed = FALSE)  # node betweenness: number of shortest paths going through a node
nodebetweenness <- sort(nodebetweenness, decreasing = TRUE) # same when values are included in g_diff or not
top100genes <- names(nodebetweenness[1:100])

# Plot network properties in one plot
data.frame("nodedegree" = nodedegree,
           "nodebetweenness" = nodebetweenness) %>%
  tidyr::gather(key = "property", value = "value", nodedegree:nodebetweenness) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~property, scales = "free") +
  xlab("") +
  ggtitle(paste0("Baseline expression network in " ,region, " (",
                 gorder(g_base), " nodes, ", gsize(g_base), " edges)"))
ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                         region,"_funcoup_focusHubGenes_baselineNetwork_betacutoff",beta_cutoff,".png"),
       width = 8, height = 6)


# # 6. Compare top genes in our network to top genes in prior
# # 6a. Rank plot to compare ranks of top genes according to nodebetweenness
# png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                       region,"_funcoup_focusHubGenes_rankPlot_nodebetweenness_prior-betacutoff",beta_cutoff,".png"),
#     width = 700, height = 700)
# plotRanks(top100genes_prior, top100genes,
#           labels = FALSE)
# dev.off()
# # Rank plot to compare ranks of top genes according to nodebetweenness
# png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                       region,"_funcoup_focusHubGenes_rankPlot_nodedegree_prior-betacutoff",beta_cutoff,".png"),
#     width = 700, height = 700)
# plotRanks(names(nodedegree_prior)[1:100], names(nodedegree)[1:100],
#           labels = FALSE)
# dev.off()
# 
# # 6b. Venn diagram (euler plot) to compare overlap of top genes according to nodebetweenness
# list1 <- list(top100genes_prior, 
#               top100genes)
# names(list1) <- c("Funcoup prior",
#                   paste("Beta cutoff", beta_cutoff))
# png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                       region,"_funcoup_focusHubGenes_vennEuler_nodebetweenness_prior-betacutoff",beta_cutoff,".png"),
#     width = 700, height = 700)
# print(plot(euler(list1, shape = "ellipse"), 
#            labels = list(cex = 1.5), quantities = list(cex = 1.5)))
# dev.off()
# # Venn diagram (euler plot) to compare overlap of top genes according to nodedegree
# list1 <- list(names(nodedegree_prior)[1:100], 
#               names(nodedegree)[1:100])
# names(list1) <- c("Funcoup prior",
#                   paste("Beta cutoff", beta_cutoff))
# png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                       region,"_funcoup_focusHubGenes_vennEuler_nodedegree_prior-betacutoff",beta_cutoff,".png"),
#     width = 700, height = 700)
# print(plot(euler(list1, shape = "ellipse"), 
#            labels = list(cex = 1.5), quantities = list(cex = 1.5)))
# dev.off()


# 6a. Compare estimated baseline network to prior (ranks of nodedegree/betweenness)
# rank of genes in prior
nodedegree_prior_rank <- rank(-nodedegree_prior[names(nodedegree)])
nodebetweenness_prior_rank <- rank(-nodebetweenness_prior[names(nodebetweenness)])
# rank of genes in diff network
nodedegree_rank <- rank(-nodedegree)
nodebetweenness_rank <- rank(-nodebetweenness)

# Spearman's rank correlation
cor_nodedegree <- cor.test(nodedegree_prior_rank, nodedegree_rank,
                           method = "spearman")
cor_nodebetweenness <- cor.test(nodebetweenness_prior_rank, nodebetweenness_rank,
                                method = "spearman")

# Scatterplot between prior and base network with correlation in label
# data.frame("prior" = nodedegree_prior[names(nodedegree)],
#            "differential" = nodedegree) %>%
#   ggplot(aes(x=prior, y=differential)) +
#   # geom_point(size=1,alpha = 0.1)
#   geom_hex() +
#   # geom_bin2d()
#   xlab("nodedegree prior network") +
#   ylab("nodedegree differential network") +
#   ggtitle(paste0("Nodedegree in prior and differential network (Correlation between ranks: ", 
#                  round(cor_nodedegree$estimate, digits = 2),")"))
# ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                          region,"_prior_differential_funcoup_correlationNodedegree_betacutoff",beta_cutoff,".png"),
#        width = 9, height = 8)
# 
# data.frame("prior" = nodebetweenness_prior[names(nodebetweenness)],
#            "differential" = nodebetweenness) %>%
#   ggplot(aes(x=prior, y=differential)) +
#   # geom_point(size=1,alpha = 0.1)
#   geom_hex() +
#   # geom_bin2d()
#   xlab("nodebetweenness prior network") +
#   ylab("nodebetweenness differential network") +
#   ggtitle(paste0("Nodebetweenness in prior and differential network (Correlation between ranks: ", 
#                  round(cor_nodebetweenness$estimate, digits = 2),")"))
# ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                          region,"_prior_differential_funcoup_correlationNodebetweenness_betacutoff",beta_cutoff,".png"),
#        width = 9, height = 8)



# tmp <- relations_prior[relations_prior$from == "ENSMUSG00000021660" | relations_prior$to == "ENSMUSG00000021660",]
# tmp1 <- relations[relations$from == "ENSMUSG00000021660" | relations$to == "ENSMUSG00000021660",]
# tmp_join <- inner_join(tmp, tmp1, 
#                    by = c("from", "to"),
#                    suffix = c(".prior", ".diff"))
# dex_notBase <- anti_join(data_dex1, data_dex0, by = c("target", "predictor"))


# # 8. Subset igraph object to top genes with their neighbours
# g1 <- induced.subgraph(graph=g_diff,vids=unlist(neighborhood(graph=g_diff,order=1,nodes=c(top100genes[1]))))
# png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                        region,"_funcoup_focusHubGenes_networkTopGene1_betacutoff",beta_cutoff,".png"),
#      width = 1000, height = 1000)
# plot(g1)
# dev.off()
# library(RCy3)
# cytoscapePing()
# createNetworkFromIgraph(g1,"myIgraph")


# 9. "Normalize" nodebetweenness by nodebetweenness in prior
# check if nodebetweenness can not be compared like this because here we use z scores as value
# and in prior no values are included
nodebetweenness_norm <- nodebetweenness/nodebetweenness_prior[names(nodebetweenness)]
nodebetweenness_mat <- data.frame("nodebetweenness" = nodebetweenness, 
                                  "nodebetweenness_prior" = nodebetweenness_prior[names(nodebetweenness)], 
                                  "nodebetweenness_norm" = nodebetweenness_norm) 
# set norm nodebetweenness to NA if nodebetweenness in prior < 10000
nodebetweenness_mat$nodebetweenness_norm[nodebetweenness_mat$nodebetweenness_prior < 10000] <- NA
nodebetweenness_mat <- arrange(nodebetweenness_mat, desc(nodebetweenness_norm))
# add column with gene symbol
nodebetweenness_mat$gene_symbol <- mapIds(org.Mm.eg.db, keys = rownames(nodebetweenness_mat), 
                                          keytype = "ENSEMBL", column="SYMBOL")
# rank normalized nodebetweenness
nodebetweenness_norm_rank <- rank(-nodebetweenness_norm[names(nodebetweenness)])

# correlation between ranks of prior nodebetweenness and norm betweenness
cor_nodebetweenness_norm <- cor.test(nodebetweenness_prior_rank, nodebetweenness_norm_rank,
                                     method = "spearman")

data.frame("prior" = nodebetweenness_prior[names(nodebetweenness)],
           "baseline" = nodebetweenness_norm[(names(nodebetweenness))]) %>%
  ggplot(aes(x=prior, y=baseline)) +
  # geom_point(size=1,alpha = 0.1)
  # geom_hex() +
  geom_bin2d() +
  xlab("nodebetweenness prior network") +
  ylab("norm. nodebetweenness baseline network") +
  ggtitle(paste0("Nodebetweenness in prior and baseline network (Correlation between ranks: ", 
                 round(cor_nodebetweenness_norm$estimate, digits = 2),")"))
ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                         region,"_prior_baseline_funcoup_correlationNodebetweennessNorm_betacutoff",beta_cutoff,".png"),
       width = 9, height = 8)

# nodebetweenness_norm <- sort(nodebetweenness_norm, decreasing = TRUE)
# top100genes_between_norm <- names(nodebetweenness_norm)[1:100]

# write nodebetweenness table to file
nodebetweenness_mat <- tibble::rownames_to_column(nodebetweenness_mat, "ensembl_id")
fwrite(nodebetweenness_mat, file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                          "03_",region,"_funcoup_baseline_nodebetweennessNorm_betacutoff",beta_cutoff,".csv"),
       quote = FALSE)


# 10. "Normalize" nodedegree by nodedegree in prior
nodedegree_norm <- nodedegree/nodedegree_prior[names(nodedegree)]
nodedegree_mat <- data.frame("nodedegree" = nodedegree, "nodedegree_prior" = nodedegree_prior[names(nodedegree)], 
                             "nodedegree_norm" = nodedegree_norm) %>%
  arrange(desc(nodedegree_norm))
# set norm nodebetweenness to NA if nodedegree in prior < 50
nodedegree_mat$nodedegree_norm[nodedegree_mat$nodedegree_prior < 50] <- NA
nodedegree_mat <- arrange(nodedegree_mat, desc(nodedegree_norm))
# add column with gene symbol
nodedegree_mat$gene_symbol <- mapIds(org.Mm.eg.db, keys = rownames(nodedegree_mat), 
                                     keytype = "ENSEMBL", column="SYMBOL")

# rank normalized nodedegree
nodedegree_norm_rank <- rank(-nodedegree_norm[names(nodedegree)])

# correlation between ranks of prior nodebetweenness and norm betweenness
cor_nodedegree_norm <- cor.test(nodedegree_prior_rank, nodedegree_norm_rank,
                                method = "spearman")

data.frame("prior" = nodedegree_prior[names(nodedegree)],
           "baseline" = nodedegree_norm[(names(nodedegree))]) %>%
  ggplot(aes(x=prior, y=baseline)) +
  # geom_point(size=1,alpha = 0.1)
  # geom_hex() +
  geom_bin2d() +
  xlab("nodedegree prior network") +
  ylab("norm. nodedegree baseline network") +
  ggtitle(paste0("Nodedegree in prior and baseline network (Correlation between ranks: ", 
                 round(cor_nodedegree_norm$estimate, digits = 2),")"))
ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                         region,"_prior_baseline_funcoup_correlationNodedegreeNorm_betacutoff",beta_cutoff,".png"),
       width = 9, height = 8)

# write nodedegree table to file
nodedegree_mat <- tibble::rownames_to_column(nodedegree_mat, "ensembl_id")
fwrite(nodedegree_mat, file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                     "03_",region,"_funcoup_baseline_nodedegreesNorm_betacutoff",beta_cutoff,".csv"),
       quote = FALSE)


# # checkout what happens to FKBP5 (ENSMUSG00000024222) in network 
# g2 <- induced.subgraph(graph=g_diff,vids=unlist(neighborhood(graph=g_diff,order=1,nodes=c("ENSMUSG00000024222"))))

