##################################################
## Project: DexStim Mouse Brain
## Date: 26.10.2020
## Author: Nathalie
##################################################
# Analyze partial correlations and create network
# Separate on dex and baseline

library(data.table)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)
library(igraph)
library(splineTimeR)

reg <- "PVN"
d <- 0

basepath <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
parcor <- paste0(basepath,"tables/coExpression_kimono/04_parCorr_singleRegion_",reg,"_dex",d,".csv")
figure_path <- paste0(basepath, "figures/02_CoExp_Kimono/01_CutoffSelection/")
GR_genes <- paste0(basepath,"data/kimono_input/63genes_ZimmermannPaper.csv")
biogrid_file <- paste0(basepath, "data/kimono_input/prior_expr_biogrid_mm.csv")

# 1. Read data
data <- fread(parcor)
hist(data$pcor)
min(abs(data$pcor))


# 2. Make network statistics for different correlation cutoffs
cutoff_vec <- seq(from = 0, to = 0.0015, by = 0.0001)
nr_edges <- rep(x = 0, length(cutoff_vec))
nr_nodes <- rep(x = 0, length(cutoff_vec))
degrees <- data.frame(cutoff = numeric(),
                      node = character(),
                      degree = numeric())
link_density <- rep(x = 0, length(cutoff_vec))
for(i in 1:length(cutoff_vec)){
  
  subset <- data[abs(data$pcor) >= cutoff_vec[i],]
  # create an igraph network from dataframe
  actors<-unique(c(subset$node1_name,subset$node2_name))
  relations <- data.frame(from=subset$node1_name,
                          to=subset$node2_name,
                          value=subset$pcor,
                          performance=subset$pval) 
  g_subset <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  nr_edges[i] <- gsize(g_subset)  # number of edges in network
  nr_nodes[i] <- gorder(g_subset) # number of nodes in network
  link_density[i] <- edge_density(g_subset)   # ratio of the number of edges and the number of possible edges
  nodedegree <- igraph::degree(g_subset)     # node degree for each node in network
  degrees <- rbind(degrees, data.frame(cutoff = rep(cutoff_vec[i], length(nodedegree)), node = names(nodedegree), degree = nodedegree))
  
}

names(nr_edges) <- cutoff_vec
png(filename = paste0(figure_path,"parCorr_nr_edges_",reg,".png"), width = 800, height = 600)
barplot(nr_edges, 
        xlab = "cutoff",
        ylab = "number of edges",
        main = paste0("Number of edges in partial correlation network: ", reg))
dev.off()
names(nr_nodes) <- cutoff_vec
png(filename = paste0(figure_path,"parCorr_nr_nodes_",reg,".png"), width = 800, height = 600)
barplot(nr_nodes,
        xlab = "cutoff",
        ylab = "number of nodes",
        main = paste0("Number of nodes in partial correlation network: ", reg))
dev.off()
names(link_density) <- cutoff_vec
png(filename = paste0(figure_path,"parCorr_link_density_",reg,".png"), width = 800, height = 600)
barplot(link_density,
        xlab = "cutoff",
        ylab = "link density",
        main = paste0("Link density in partial correlation network: ", reg))
dev.off()

ggplot(degrees, aes(x = degree)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 100) +
  facet_wrap(. ~ cutoff, ncol = 5) +
  ggtitle(paste0("Distribution of node degrees for different cutoffs in partial correlation network: ", reg))
ggsave(filename = paste0(figure_path,"parCorr_node_degrees_",reg,".png"), width = 10, height = 7)


# 3. Calculate betweenness and fit powerlaw for different correlation cutoffs
cutoff_vec <- seq(from = 0.0005, to = 0.0008, by = 0.0001)
between <- data.frame(cutoff = numeric(),
                          node = character(),
                          betweenness = numeric())
powerlaw_exponent <- rep(x = 0, length(cutoff_vec))
for(i in 1:length(cutoff_vec)){
  
  subset <- data[abs(data$pcor) >= cutoff_vec[i],]
  # create an igraph network from dataframe
  actors<-unique(c(subset$node1_name,subset$node2_name))
  relations <- data.frame(from=subset$node1_name,
                          to=subset$node2_name,
                          value=subset$pcor,
                          performance=subset$pval) 
  g_subset <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  nodebetweenness <- betweenness(g_subset, directed = FALSE)  # node betweenness: number of shortest paths going through a node
  between <- rbind(between, data.frame(cutoff = rep(cutoff_vec[i], length(nodebetweenness)), node = names(nodebetweenness), betweenness = nodebetweenness))
  
  # fit a scale free network to data
  scaleFreeProp <- networkProperties(g_subset)
  powerlaw_exponent[i] <- scaleFreeProp[1,3]
  file.rename(from="scale_free_properties.pdf", to = paste0(figure_path,"parCorr_scale_free_properties_",reg,"_",cutoff_vec[i],".pdf"))
}

ggplot(between, aes(x = betweenness)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 100) +
  facet_wrap(. ~ cutoff, ncol = 5) +
  ggtitle(paste0("Distribution of node betweenness for different cutoffs in partial correlation network: ", reg))
ggsave(filename = paste0(figure_path,"parCorr_node_betweenness_",reg,".png"), width = 10, height = 7)

names(powerlaw_exponent) <- cutoff_vec
png(filename = paste0(figure_path,"parCorr_powerlaw_exponent_",reg,".png"), width = 800, height = 600)
barplot(powerlaw_exponent,
        xlab = "cutoff",
        ylab = "powerlaw_exponent",
        main = paste0("Power-law exponent in partial correlation network: ", reg))
dev.off()



# # 4. Read BioGrid as reference
# biogrid <- fread(biogrid_file)
# 
# cutoff_vec <- seq(from = 0, to = 0.0015, by = 0.0001)
# fisher_pval <- rep(0, length(cutoff_vec))
# for(i in 1:length(cutoff_vec)){
#   subset <- data[abs(data$pcor) >= cutoff_vec[i],]
#   
#   # Identify nodes present in Biogrid and parCorr subnet
#   nodes_subset <- unique(c(subset$node1_name,subset$node2_name))
#   nodes_biogrid <- unique(c(biogrid$ensembl_A, biogrid$ensembl_B))
#   inters_biosub <- intersect(nodes_subset, nodes_biogrid)
#   
#   # Subset biogrid edges and subset edges to these nodes
#   subset_biogrid <- biogrid %>%
#     filter(ensembl_A %in% inters_biosub & ensembl_B %in% inters_biosub)
#   subset_parCorr <- subset %>%
#     filter(node1_name %in% inters_biosub & node2_name %in% inters_biosub) %>%
#     mutate(from=node1_name, to=node2_name)
#   
#   # Create adjacency matrix for biogrid subnet und subset subnet
#   g_biograph <- graph_from_data_frame(subset_biogrid %>% mutate(from=ensembl_A, to=ensembl_B),
#                                       directed=FALSE, vertices = inters_biosub)
#   a_biograph <- as_adjacency_matrix(g_biograph,names=TRUE,sparse=FALSE,type='lower')
#   
#   g_subset <- graph_from_data_frame(data.frame(subset_parCorr$from, subset_parCorr$to),
#                                     directed = FALSE, vertices = inters_biosub)
#   a_subset <- as_adjacency_matrix(g_subset, names = TRUE, sparse = FALSE, type = 'lower')
#   
#   biogrid_vector <- a_biograph[lower.tri(a_biograph)]
#   subset_vector <- a_subset[lower.tri(a_subset)]
#   
#   contingency <- table(biogrid_vector, subset_vector)
#   f <- fisher.test(contingency)
#   fisher_pval[i] <- f$p.value
#   
# }



# Read 63 genes from Zimmermann/Arloth Paper
GRgenes <- fread(GR_genes)

cutoff_vec <- seq(from = 0, to = 0.0015, by = 0.0001)
nr_edges <- rep(x = 0, length(cutoff_vec))
nr_nodes <- rep(x = 0, length(cutoff_vec))
for(i in 1:length(cutoff_vec)){
  
  subset <- data[abs(data$pcor) >= cutoff_vec[i],]
  subset <- subset[node1_name %in% GRgenes$Ensembl,]
  subset <- subset[node2_name %in% GRgenes$Ensembl , ]
  # create an igraph network from dataframe
  actors<-unique(c(subset$node1_name,subset$node2_name))
  relations <- data.frame(from=subset$node1_name,
                          to=subset$node2_name,
                          value=subset$pcor,
                          performance=subset$pval) 
  g_subset <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  nr_edges[i] <- gsize(g_subset)  # number of edges in network
  nr_nodes[i] <- gorder(g_subset) # number of nodes in network
  
}

names(nr_edges) <- cutoff_vec
png(filename = paste0(figure_path,"parCorr_63genes_nr_edges_",reg,".png"), width = 800, height = 600)
barplot(nr_edges, 
        xlab = "cutoff",
        ylab = "number of edges",
        main = paste0("Number of edges in partial correlation network (63 genes): ", reg))
dev.off()
names(nr_nodes) <- cutoff_vec
png(filename = paste0(figure_path,"parCorr_63genes_nr_nodes_",reg,".png"), width = 800, height = 600)
barplot(nr_nodes,
        xlab = "cutoff",
        ylab = "number of nodes",
        main = paste0("Number of nodes in partial correlation network (63 genes): ", reg))
dev.off()
