##################################################
## Project: DexStim Mouse Brain
## Date: 26.10.2020
## Author: Nathalie
##################################################
# Analyze kimono output and create network
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
kimono <- paste0(basepath,"tables/coExpression_kimono/04_singleRegion_",reg,"_dex",d,"_funcoup.csv")
figure_path <- paste0(basepath, "figures/02_CoExp_Kimono/01_CutoffSelection/")
GR_genes <- paste0(basepath,"data/kimono_input/63genes_ZimmermannPaper.csv")
#biogrid_file <- paste0(basepath, "data/kimono_input/prior_expr_biogrid_mm.csv")

# 1. Read data
data <- fread(kimono)
hist(data$value, breaks = 40000,
     xlim = c(-0.001,0.001))
max(abs(data$value))
hist(data$performance)
data <- data %>%
  filter(performance >= 0.1)

# 2. Make network statistics for different correlation cutoffs
cutoff_vec <- c(0,0.000001,0.00001,0.0001,0.0005,0.001,0.002,0.005,0.01,0.1)
nr_edges <- rep(x = 0, length(cutoff_vec))
nr_nodes <- rep(x = 0, length(cutoff_vec))
degrees <- data.frame(cutoff = numeric(),
                      node = character(),
                      degree = numeric())
link_density <- rep(x = 0, length(cutoff_vec))
for(i in 1:length(cutoff_vec)){
  
  subset <- data[abs(data$value) >= cutoff_vec[i],]
  # create an igraph network from dataframe
  actors<-unique(c(subset$target,subset$predictor))
  relations <- data.frame(from=subset$target,
                          to=subset$predictor,
                          value=subset$value,
                          performance=subset$performance) 
  g_subset <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  nr_edges[i] <- gsize(g_subset)  # number of edges in network
  nr_nodes[i] <- gorder(g_subset) # number of nodes in network
  link_density[i] <- edge_density(g_subset)   # ratio of the number of edges and the number of possible edges
  nodedegree <- igraph::degree(g_subset)     # node degree for each node in network
  degrees <- rbind(degrees, data.frame(cutoff = rep(cutoff_vec[i], length(nodedegree)), node = names(nodedegree), degree = nodedegree))
  
}

names(nr_edges) <- cutoff_vec
png(filename = paste0(figure_path,"kimono_nr_edges_",reg,".png"), width = 800, height = 600)
barplot(nr_edges, 
        xlab = "cutoff",
        ylab = "number of edges",
        main = paste0("Number of edges in kimono network: ", reg))
dev.off()
names(nr_nodes) <- cutoff_vec
png(filename = paste0(figure_path,"kimono_nr_nodes_",reg,".png"), width = 800, height = 600)
barplot(nr_nodes,
        xlab = "cutoff",
        ylab = "number of nodes",
        main = paste0("Number of nodes in kimono network: ", reg))
dev.off()
names(link_density) <- cutoff_vec
png(filename = paste0(figure_path,"kimono_link_density_",reg,".png"), width = 800, height = 600)
barplot(link_density,
        xlab = "cutoff",
        ylab = "link density",
        main = paste0("Link density in kimono network: ", reg))
dev.off()

ggplot(degrees, aes(x = degree)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 100) +
  facet_wrap(. ~ cutoff, ncol = 5) +
  ggtitle(paste0("Distribution of node degrees for different cutoffs in kimono network: ", reg))
ggsave(filename = paste0(figure_path,"kimono_node_degrees_",reg,".png"), width = 10, height = 7)


# 3. Calculate betweenness and fit powerlaw for different correlation cutoffs
cutoff_vec <- c(0.0001,0.0005,0.001,0.002,0.005)
between <- data.frame(cutoff = numeric(),
                      node = character(),
                      betweenness = numeric())
powerlaw_exponent <- rep(x = 0, length(cutoff_vec))
for(i in 1:length(cutoff_vec)){
  
  subset <- data[abs(data$value) >= cutoff_vec[i],]
  # create an igraph network from dataframe
  actors<-unique(c(subset$target,subset$predictor))
  relations <- data.frame(from=subset$target,
                          to=subset$predictor,
                          value=subset$value,
                          performance=subset$performance) 
  g_subset <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  nodebetweenness <- betweenness(g_subset, directed = FALSE)  # node betweenness: number of shortest paths going through a node
  between <- rbind(between, data.frame(cutoff = rep(cutoff_vec[i], length(nodebetweenness)), node = names(nodebetweenness), betweenness = nodebetweenness))
  
  # fit a scale free network to data
  scaleFreeProp <- networkProperties(g_subset)
  powerlaw_exponent[i] <- scaleFreeProp[1,3]
  file.rename(from="scale_free_properties.pdf", to = paste0(figure_path,"kimono_scale_free_properties_",reg,"_",cutoff_vec[i],".pdf"))
}

ggplot(between, aes(x = betweenness)) +
  geom_histogram(position="identity", colour="grey40", alpha=0.2, bins = 100) +
  facet_wrap(. ~ cutoff, ncol = 5) +
  ggtitle(paste0("Distribution of node betweenness for different cutoffs in kimono network: ", reg))
ggsave(filename = paste0(figure_path,"kimono_node_betweenness_",reg,".png"), width = 10, height = 7)

names(powerlaw_exponent) <- cutoff_vec
png(filename = paste0(figure_path,"kimono_powerlaw_exponent_",reg,".png"), width = 800, height = 600)
barplot(powerlaw_exponent,
        xlab = "cutoff",
        ylab = "powerlaw_exponent",
        main = paste0("Power-law exponent in kimono network: ", reg))
dev.off()



# Read 63 genes from Zimmermann/Arloth Paper
GRgenes <- fread(GR_genes)

cutoff_vec <- c(0,0.000001,0.00001,0.0001,0.0005,0.001,0.002,0.005,0.01,0.1)
nr_edges <- rep(x = 0, length(cutoff_vec))
nr_nodes <- rep(x = 0, length(cutoff_vec))
for(i in 1:length(cutoff_vec)){
  
  subset <- data[abs(data$value) >= cutoff_vec[i],]
  subset <- subset[target %in% GRgenes$Ensembl,]
  subset <- subset[predictor %in% GRgenes$Ensembl , ]
  # create an igraph network from dataframe
  actors<-unique(c(subset$target,subset$predictor))
  relations <- data.frame(from=subset$target,
                          to=subset$predictor,
                          value=subset$value,
                          performance=subset$performance) 
  g_subset <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  nr_edges[i] <- gsize(g_subset)  # number of edges in network
  nr_nodes[i] <- gorder(g_subset) # number of nodes in network
  
}

names(nr_edges) <- cutoff_vec
png(filename = paste0(figure_path,"kimono_63genes_nr_edges_",reg,".png"), width = 800, height = 600)
barplot(nr_edges, 
        xlab = "cutoff",
        ylab = "number of edges",
        main = paste0("Number of edges in kimono network (63 genes): ", reg))
dev.off()
names(nr_nodes) <- cutoff_vec
png(filename = paste0(figure_path,"kimono_63genes_nr_nodes_",reg,".png"), width = 800, height = 600)
barplot(nr_nodes,
        xlab = "cutoff",
        ylab = "number of nodes",
        main = paste0("Number of nodes in kimono network (63 genes): ", reg))
dev.off()

