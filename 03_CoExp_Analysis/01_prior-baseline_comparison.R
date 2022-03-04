##################################################
## Project: DexStim Mouse Brain
## Date: 12.12.2020
## Author: Nathalie
##################################################
# Compare prior and baseline network (nodedegree and nodebetweenness)

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(WGCNA)
library(eulerr)

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

### ANALYSIS ---------------------------------------

# 1. Read data
# 1a. Read co expression networks
dex0_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_singleRegion\\_",region,"\\_dex0\\_funcoup\\_SE\\_.*\\.csv"),
                         full.names = TRUE)
# dex1_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
#                          pattern = paste0("04\\_singleRegion\\_",region,"\\_dex1\\_funcoup\\_SE\\_.*\\.csv"),
#                          full.names = TRUE)

data_dex0 <- readFiles_concat(dex0_files)
# data_dex1 <- readFiles_concat(dex1_files)

# 1b. Read prior and make network
funcoup_prior <- fread(file = paste0(basepath, "data/kimono_input/prior_expr_funcoup_mm.csv"))
nodes_prior <- unique(c(funcoup_prior$Gene_A, funcoup_prior$Gene_B))
relations_prior <- data.frame(from = funcoup_prior$Gene_A,
                              to = funcoup_prior$Gene_B)
g_prior <- graph_from_data_frame(relations_prior, directed=FALSE, vertices=nodes_prior)
# Calculate nodedegrees of prior
nodedegree_prior <- igraph::degree(g_prior)
nodedegree_prior <- sort(nodedegree_prior, decreasing = TRUE)
# Calculate nodebetweenness of prior
# nodebetweenness_prior <- betweenness(g_prior, directed = FALSE)  # node betweenness: number of shortest paths going through a node
# saveRDS(nodebetweenness_prior, file = paste0(basepath, "data/workspaces/nodebetweenness_funcoup.rds"))
nodebetweenness_prior <- readRDS(paste0(basepath, "data/workspaces/nodebetweenness_funcoup.rds"))
nodebetweenness_prior <- sort(nodebetweenness_prior, decreasing = TRUE)

# 2. Remove interactions with very low r squared values & intercept & SVs
data <- data_dex0 %>% 
  filter(overall_rsq >= rsquared_cutoff) %>%
  filter(predictor != '(Intercept)')
data <- data[!startsWith(data$predictor, "SV"),]
# remove duplicated interactions (mistake made when separating nodes into chunks)
data <- data %>%
  distinct(target, predictor, .keep_all = TRUE)

# 3. Keep only interactions that a beta value > cutoff
data <- data %>%
  mutate(high_beta = (abs(beta_mean) > beta_cutoff))
data_cut <- data %>%
  filter(high_beta)

# 4. Create network
# Nodes in network
node_vec <- unique(c(data_cut$target, data_cut$predictor))

# Find modules
relations <- data.frame(from=data_cut$target,
                        to=data_cut$predictor,
                        value=data_cut$beta_mean)
# relations <- data.frame(from=data_diff$target,
#                         to=data_diff$predictor) 
g <- graph_from_data_frame(relations, directed=FALSE, vertices=node_vec)
# does graph contain multiple edges with same start and endpoint or loop edges
is_simple(g)
g <- simplify(g) # check the edge attribute parameter

# Calculate nodedegree
nodedegree <- igraph::degree(g)
nodedegree <- sort(nodedegree, decreasing = TRUE)
nodedegree_rank <- rank(-nodedegree)

# Calculate nodebetweenness
nodebetweenness <- betweenness(g, directed = FALSE)  # node betweenness: number of shortest paths going through a node
nodebetweenness <- sort(nodebetweenness, decreasing = TRUE) # same when values are included in g_diff or not
nodebetweenness_rank <- rank(-nodebetweenness)


# 5. Compare estimated baseline network to prior
# rank of genes in prior
nodedegree_prior_rank <- rank(-nodedegree_prior[names(nodedegree)])
nodebetweenness_prior_rank <- rank(-nodebetweenness_prior[names(nodebetweenness)])

# Spearman's rank correlation
cor_nodedegree <- cor.test(nodedegree_prior_rank, nodedegree_rank,
                           method = "spearman")
cor_nodebetweenness <- cor.test(nodebetweenness_prior_rank, nodebetweenness_rank,
                                method = "spearman")

# Scatterplot between prior and base network with correlation in label
data.frame("prior" = nodedegree_prior[names(nodedegree)],
           "baseline" = nodedegree) %>%
ggplot(aes(x=prior, y=baseline)) +
  # geom_point(size=1,alpha = 0.1)
  geom_hex() +
  # geom_bin2d()
  xlab("nodedegree prior network") +
  ylab("nodedegree baseline network") +
  ggtitle(paste0("Nodedegree in prior and baseline network (Correlation between ranks: ", 
                 round(cor_nodedegree$estimate, digits = 2),")"))
ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                         region,"_prior_baseline_funcoup_correlationNodedegree_betacutoff",beta_cutoff,".png"),
       width = 9, height = 8)

data.frame("prior" = nodebetweenness_prior[names(nodebetweenness)],
           "baseline" = nodebetweenness) %>%
  ggplot(aes(x=prior, y=baseline)) +
  # geom_point(size=1,alpha = 0.1)
  geom_hex() +
  # geom_bin2d()
  xlab("nodebetweenness prior network") +
  ylab("nodebetweenness baseline network") +
  ggtitle(paste0("Nodebetweenness in prior and baseline network (Correlation between ranks: ", 
                 round(cor_nodebetweenness$estimate, digits = 2),")"))
ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                         region,"_prior_baseline_funcoup_correlationNodebetweenness_betacutoff",beta_cutoff,".png"),
       width = 9, height = 8)






