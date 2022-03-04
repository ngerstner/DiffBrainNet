##################################################
## Project: DexStim Mouse Brain
## Date: 21.01.2020
## Author: Nathalie
##################################################
# Analyze multitissue network
# for comparison with single tissue networks

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(eulerr)
library(org.Mm.eg.db)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
regions <- c("AMY", "CER", "dCA1", "dDG", "PFC", "PVN", "vCA1", "vDG")
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

# subset network to only one brain region
subset_network <- function(network, region){
  
  # find regions that have connections to brain region
  select_genes <- network$target[network$predictor == region]
  
  # subset network data to those targets with connection to region
  network <- network %>%
    filter(target %in% select_genes) %>%
    filter(startsWith(predictor, "ENSMUS"))
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


# 2a. Read kimono expression networks
dex0_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_multipleRegions\\_funcoup\\_dex0\\_SE\\_.*\\.csv"),
                         full.names = TRUE)
dex1_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_multipleRegions\\_funcoup\\_dex1\\_SE\\_.*\\.csv"),
                         full.names = TRUE)

data_dex0 <- readFiles_concat(dex0_files)
data_dex1 <- readFiles_concat(dex1_files)

# 2b. Join Base and Dex data frame
data <- inner_join(data_dex0, data_dex1, 
                   by = c("target", "predictor"),
                   suffix = c(".base", ".dex"))
dex_notBase <- anti_join(data_dex1, data_dex0, by = c("target", "predictor"))
head(data)
rm(data_dex0)
rm(data_dex1)


# 3. Remove interactions with very low r squared values & intercept & SVs
data <- data %>% 
  filter(overall_rsq.base >= rsquared_cutoff, overall_rsq.dex >= rsquared_cutoff) %>%
  filter(predictor != '(Intercept)')
data <- data[!startsWith(data$predictor, "SV"),]
# remove duplicated interactions (mistake made when separating nodes into chunks)
data <- data %>%
  distinct(target, predictor, .keep_all = TRUE)


# 4a. Calculate z scores for interactions that are left
data <- mutate(data, z = z_score(beta_mean.dex,beta_mean.base, beta_stderr.dex, beta_stderr.base))
hist(data$beta_mean.base)

# 4b. Keep only interactions that have at least in one network a beta value > cutoff
data <- data %>%
  mutate(diff = (abs(beta_mean.base) > beta_cutoff | abs(beta_mean.dex) > beta_cutoff))
data_diff1 <- data %>%
  filter(diff)
# calculate p-value for z-score
data_diff1$p_diff <- 2*pnorm(-abs(data_diff1$z))
data_diff1$p_adj <- p.adjust(data_diff1$p_diff, method = "fdr")


# 5. Create diff network corresponding to beta cutoff
data_diff <- data_diff1 %>%
  filter(p_adj <= padj_cutoff)
head(data_diff[,c("beta_mean.base", "beta_stderr.base", "beta_mean.dex", "beta_stderr.dex", "z", "p_adj")], 20)

# Nodes in network
node_vec <- unique(c(data_diff$target, data_diff$predictor))
node_type <- startsWith(node_vec, "ENSMUS")
names(node_type) <- node_vec

# Edges in network
relations <- data.frame(from=data_diff$target,
                        to=data_diff$predictor,
                        value=data_diff$z,
                        performance=data_diff$p_adj)
# relations <- data.frame(from=data_diff$target,
#                         to=data_diff$predictor) 
g_diff <- graph_from_data_frame(relations, directed=FALSE, vertices=node_vec)
# does graph contain multiple edges with same start and endpoint or loop edges
is_simple(g_diff)
g_diff <- simplify(g_diff) # check the edge attribute parameter


# Calculate nodedegree
nodedegree <- igraph::degree(g_diff)
nodedegree <- sort(nodedegree, decreasing = TRUE)

# Calculate nodebetweenness
nodebetweenness <- betweenness(g_diff, directed = FALSE)  # node betweenness: number of shortest paths going through a node
nodebetweenness <- sort(nodebetweenness, decreasing = TRUE) # same when values are included in g_diff or not
top100genes <- names(nodebetweenness[1:100])

# Plot network properties in one plot
data.frame("nodedegree" = nodedegree,
           "nodebetweenness" = nodebetweenness,
           "gene" = node_type[names(nodedegree)]) %>%
  tidyr::gather(key = "property", value = "value", nodedegree:nodebetweenness) %>%
  ggplot(aes(y = value, fill = gene)) +
  geom_boxplot() +
  facet_wrap(~property, scales = "free") +
  theme_light() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10), strip.text.x = element_text(size = 12)) +
  scale_fill_discrete(labels = c("tissues", "genes")) +
  ggtitle(paste0("Differential expression network for all brain regions (",
                 gorder(g_diff), " nodes, ", gsize(g_diff), " edges)"))
ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/multipleRegions/",
                         "08_funcoup_focusHubGenes_diffNetwork_betacutoff",beta_cutoff,".png"),
       width = 8, height = 6)


# ANALYZE EACH BRAIN REGION ------------------------

for (region in regions){
  
  # A1. Get single region network with top genes
  # Subset network
  net_region <- subset_network(data_diff, region)
  
  # Nodes in network
  node_vec <- unique(c(net_region$target, net_region$predictor))
  
  # Edges in network
  relations <- data.frame(from=net_region$target,
                          to=net_region$predictor,
                          value=net_region$z,
                          performance=net_region$p_adj)
  g_diff_reg <- graph_from_data_frame(relations, directed=FALSE, vertices=node_vec)
  # does graph contain multiple edges with same start and endpoint or loop edges
  is_simple(g_diff_reg)
  g_diff <- simplify(g_diff_reg) # check the edge attribute parameter
  
  # Calculate nodedegree
  nodedegree_reg <- igraph::degree(g_diff_reg)
  nodedegree_reg <- sort(nodedegree_reg, decreasing = TRUE)
  
  # Calculate nodebetweenness
  nodebetweenness_reg <- betweenness(g_diff_reg, directed = FALSE)  # node betweenness: number of shortest paths going through a node
  nodebetweenness_reg <- sort(nodebetweenness_reg, decreasing = TRUE) # same when values are included in g_diff or not
  top100genes_reg <- names(nodebetweenness_reg[1:100])
  
  
  # A2. Compare extracted single region network with kimono single region network
  # Read nodedegree and nodebetweenness 
  nodedegree_single <- fread(file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                        "04_",region,"_funcoup_differential_nodedegreesNorm_betacutoff",beta_cutoff,".csv"),
                          quote = FALSE)
  nodebetweenness_single <- fread(file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                           "04_",region,"_funcoup_differential_nodebetweennessNorm_betacutoff",beta_cutoff,".csv"),
                             quote = FALSE)
  
  # A3. Comparison nodedegree
  # Rank nodedegree
  nodedegree_multi_rank <- rank(-nodedegree_reg[nodedegree_single$ensembl_id])
  nodedegree_single_rank <- rank(-nodedegree_single$nodedegree)
  names(nodedegree_single_rank) <- nodedegree_single$ensembl_id

  # Correlation between ranks of nodedegrees
  cor_nodedegree <- cor.test(nodedegree_single_rank, nodedegree_multi_rank,
                                  method = "spearman")
  
  data.frame("single" = nodedegree_single$nodedegree,
             "multi" = nodedegree_reg[nodedegree_single$ensembl_id]) %>%
    ggplot(aes(x=single, y=multi)) +
    # geom_point(size=1,alpha = 0.1)
    # geom_hex() +
    geom_bin2d() +
    xlab("nodedegree single tissue network") +
    ylab("nodedegree multi tissue network") +
    ggtitle(paste0("Nodedegree in single and multi tissue network (Correlation between ranks: ", 
                   round(cor_nodedegree$estimate, digits = 2),")"))
  ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/multipleRegions/",
                           "08_",region,"_single_multi_funcoup_correlationNodedegree_betacutoff",beta_cutoff,".png"),
         width = 9, height = 8)
  
  
  # A4. Comparison nodebetweenness
  # Rank nodebetweenness
  nodebetween_multi_rank <- rank(-nodebetweenness_reg[nodebetweenness_single$ensembl_id])
  nodebetween_single_rank <- rank(-nodebetweenness_single$nodebetweenness)
  names(nodebetween_single_rank) <- nodebetweenness_single$ensembl_id
  
  # Correlation between ranks of nodedegrees
  cor_nodebetween <- cor.test(nodebetween_single_rank, nodebetween_multi_rank,
                             method = "spearman")
  
  data.frame("single" = nodebetweenness_single$nodebetweenness,
             "multi" = nodebetweenness_reg[nodebetweenness_single$ensembl_id]) %>%
    ggplot(aes(x=single, y=multi)) +
    # geom_point(size=1,alpha = 0.1)
    # geom_hex() +
    geom_bin2d() +
    xlab("nodebetweenness single tissue network") +
    ylab("nodebetweenness multi tissue network") +
    ggtitle(paste0("Nodebetweenness in single and multi tissue network (Correlation between ranks: ", 
                   round(cor_nodedegree$estimate, digits = 2),")"))
  ggsave(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/multipleRegions/",
                           "08_",region,"_single_multi_funcoup_correlationNodebetweenness_betacutoff",beta_cutoff,".png"),
         width = 9, height = 8)

}

