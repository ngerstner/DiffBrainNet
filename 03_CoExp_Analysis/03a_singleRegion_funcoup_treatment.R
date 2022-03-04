##################################################
## Project: DexStim Mouse Brain
## Date: 04.10.2021
## Author: Nathalie
##################################################
# Save filtered treatment network
# Use beta cutoff and rsquared cutoff as for baseline and differential network

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(eulerr)
library(org.Mm.eg.db)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
region <- "AMY"
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


### ANALYSIS ---------------------------------------

# 1. Read kimono treatment expression networks
dex1_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_singleRegion\\_",region,"\\_dex1\\_funcoup\\_SE\\_.*\\.csv"),
                         full.names = TRUE)

data_dex1 <- readFiles_concat(dex1_files)


# 2. Remove interactions with very low r squared values & intercept & SVs
data <- data_dex1 %>% 
  filter(overall_rsq >= rsquared_cutoff) %>%
  filter(predictor != '(Intercept)')
data <- data[!startsWith(data$predictor, "SV"),]
# remove duplicated interactions (mistake made when separating nodes into chunks)
data <- data %>%
  distinct(target, predictor, .keep_all = TRUE)


# 3. Keep only interactions that have a beta value > cutoff
data <- data %>%
  mutate(betacut = (abs(beta_mean) > beta_cutoff))
data_cut <- data %>%
  filter(betacut)


# 4. Create treatment network corresponding to beta cutoff
head(data_cut[,c("target", "predictor", "beta_mean", "beta_stderr")], 20)
# Save filtered network
fwrite(data_cut, file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                "04_singleRegion_",region,"_filtered_treatmentNetwork.csv"))


# Nodes in network
node_vec <- unique(c(data_cut$target, data_cut$predictor))

# Network properties
relations <- data.frame(from=data_cut$target,
                        to=data_cut$predictor)
g_treat <- graph_from_data_frame(relations, directed=FALSE, vertices=node_vec)
# does graph contain multiple edges with same start and endpoint or loop edges
is_simple(g_treat)
g_treat <- simplify(g_treat) # check the edge attribute parameter

# Calculate nodedegree
nodedegree <- igraph::degree(g_treat)
nodedegree <- sort(nodedegree, decreasing = TRUE)

# Calculate nodebetweenness
nodebetweenness <- betweenness(g_treat, directed = FALSE)  # node betweenness: number of shortest paths going through a node
nodebetweenness <- sort(nodebetweenness, decreasing = TRUE) # same when values are included in g_diff or not


# 1. Prior network properties for normalization
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


# 9. "Normalize" nodebetweenness by nodebetweenness in prior
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

# write nodebetweenness table to file
nodebetweenness_mat <- tibble::rownames_to_column(nodebetweenness_mat, "ensembl_id")
fwrite(nodebetweenness_mat, file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                          "03a_",region,"_funcoup_treatment_nodebetweennessNorm_betacutoff",beta_cutoff,".csv"),
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

# write nodedegree table to file
nodedegree_mat <- tibble::rownames_to_column(nodedegree_mat, "ensembl_id")
fwrite(nodedegree_mat, file = paste0(basepath, "/tables/coExpression_kimono/03_AnalysisFuncoup/",
                                     "03a_",region,"_funcoup_treatment_nodedegreesNorm_betacutoff",beta_cutoff,".csv"),
       quote = FALSE)


