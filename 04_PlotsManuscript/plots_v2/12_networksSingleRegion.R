##################################################
## Project: DexStim Mouse Brain
## Date: 29.03.2022
## Author: Nathalie
##################################################
# Abcd1 network in PFC with igraph to pdf


library(igraph)
library(data.table)
library(dplyr)
library(org.Mm.eg.db)
library(RColorBrewer)


### FUNCTIONS ------------------------------------

# Load data from one brain region
network_df <- function(network_region, network_type = "diff"){
  file_path <- paste0(
    "tables/coExpression_kimono/03_AnalysisFuncoup/",
    "04_singleRegion_",
    network_region,
    "_filtered_",
    network_type,
    "Network.csv"
  )
  # Read data
  if (network_type == "diff") {
    data <-
      fread(file = file_path) %>%
      dplyr::select(target, predictor, beta_mean.base, beta_mean.dex, z, p_adj)
  } else {
    data <-
      fread(file = file_path) %>%
      dplyr::select(target, predictor, beta_mean)
  }
  
  return(data)
}

# DE genes of region
de_genes <- function(network_region){
  # Read and subset DE genes for coloring
  df_de <- fread(paste0( "tables/02_",
                         network_region,"_deseq2_Dex_1_vs_0_lfcShrink.txt"))
  na_indices <- which(is.na(df_de$padj))
  df_de$padj[na_indices] <- 1
  df_de <- df_de[df_de$padj <= 0.1,]
  
  return(df_de$Ensembl_ID)
}

# hub genes of region
hub_genes <- function(network_region){
  # Read and subset hub genes for coloring
  df_hub <- fread(paste0("tables/coExpression_kimono/03_AnalysisFuncoup/04_",
                         network_region,"_funcoup_differential",
                         "_nodebetweennessNorm_betacutoff0.01.csv")) %>%
    filter(nodebetweenness_norm >= 1)
  
  return(df_hub$ensembl_id)
}

# bring input genes to correct format
reformat_genes <- function(list_genes){
  if (!startsWith(list_genes[1], "ENSMUSG")){
    format_gene <- sapply(list_genes, stringr::str_to_title)
    format_gene <- mapIds(org.Mm.eg.db, keys = format_gene, column = "ENSEMBL", keytype = "SYMBOL")
  } else {
    format_gene <- list_genes
  }
  return(format_gene)
}

# subset network to gene of interest and if required also
# neighbours of genes
subset_network <- function(network_dt, subset_gene, neighbours){
  
  if (neighbours) {
    e_interest <- network_dt %>% 
      filter(from %in% subset_gene | to %in% subset_gene)
    e_interest <- unique(c(e_interest$from, e_interest$to))
    network_subset <- network_dt %>%
      filter(from %in% e_interest & to %in% e_interest)
  } else {
    network_subset <- network_dt %>%
      filter(from %in% subset_gene & to %in% subset_gene)
  }
  
  return(network_subset)
}

# remove duplicated edges from network dataframe
simplify_network_df <- function(network_dt){
  
  network_dt <- network_dt %>%
    dplyr::filter(from != to) 
}

# function to generate baseline/diff network, single regions
read_network <- function(data, de_genes, hub_genes, gene, mode = "diff", 
                         id_type = "Gene Symbol", neighbours = TRUE){
  
  # input genes
  print(gene)
  
  # Edges
  if(mode == "baseline" | mode == "treatment"){
    c <- rep("grey", nrow(data))
    relations <- data.table(from=data$predictor,
                            to=data$target,
                            beta = data$beta_mean,
                            color = c)
    # width = rescale(abs(data$beta_mean), to = c(0,1)))
    ledges <- data.frame(color = c("grey"),
                         label = c("beta"), 
                         #arrows = c("right", "right"),
                         font.align = "top")
  } else {
    c <- rep("#db9512", nrow(data))
    c[data$z > 0] <- "#128bdb"
    relations <- data.table(from=data$predictor,
                            to=data$target,
                            z = data$z,
                            p_adj = data$p_adj,
                            color = c)
    ledges <- data.frame(color = c("#db9512", "#128bdb"),
                         label = c("z < 0", "z > 0"), 
                         arrows = c("right", "right"),
                         font.align = "top")
  }
  relations <- simplify_network_df(relations)
  relations <- subset_network(relations, gene, neighbours)
  print(head(relations))
  
  # Nodes
  # Generate dataframe for node layout
  node_vec <- unique(c(relations$from, relations$to))
  if (id_type == "Gene Symbol"){
    labels <- mapIds(org.Mm.eg.db, keys = node_vec,
                     column = "SYMBOL", keytype = "ENSEMBL")
  } else {
    labels <- node_vec
  }
  print(node_vec)
  
  # set colours according to DE status of gene
  col <- rep("darkblue", length(node_vec))
  col[node_vec %in% de_genes] <- "orange"
  col[node_vec %in% hub_genes] <- "darkred"
  col[node_vec %in% de_genes & node_vec %in% hub_genes] <- "#800080" # purple
  nodes <- data.table(id = node_vec,
                      label = labels,
                      color = col)
  
  # nodes data.frame for legend
  # lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
  #                      font.color = c("black", "white", "white", "white"),
  #                      shape = c( "ellipse"), 
  #                      color = c("orange", "darkred", "purple", "darkblue"),
  #                      title = "Informations", id = 1:4)
  
  # return(list("edges" = relations, "nodes" = nodes, "ledges" = ledges, "lnodes" = lnodes))
  return(list("edges" = relations, "nodes" = nodes, "ledges" = ledges))
  
}


### CREATE NETWORK Abcd1 ####################################

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"

nd <- network_df("PFC")
de <- de_genes("PFC")
hub <- hub_genes("PFC")

genes <- reformat_genes(c("Abcd1"))

n <- read_network(nd, de, hub, genes)

g <- graph_from_data_frame(n$edges, directed = TRUE, vertices = n$nodes)

# legend df for nodes
lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                     font.color = c("black", "white", "white", "white"),
                     shape = c( "ellipse"), 
                     color = c("orange", "darkred", "purple", "darkblue"),
                     title = "Informations", id = 1:4)

l <- layout_nicely(g)

pdf(file = paste0(basepath, "scripts/07_PlotsManuscript/plots_v2/12_Figure3_networkAbcd1.pdf"), 
    width = 10, height = 10)
plot(g, edge.arrow.size=.4, layout=l,
     vertex.size = 10,
     #vertex.label = ifelse(degree(g) > 4, V(g)$label, NA),
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.8,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1.1,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
     edge.arrow.size = 0.8
     )
legend(x=-1,y=-1,legend=c("DE", "hub", "no DE & no hub"),
       pch=21,
       col=c("orange", "darkred", "darkblue"), 
       pt.bg=c("orange", "darkred", "darkblue"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


### CREATE NETWORK Abc transporters pathway ####################################

# this list is not the complete Abc transporters pathway but taken over from the screenshot
genes <- reformat_genes(c("Abcd1", "Abcd3", "Pex3", "Abcb9", "Abcb4", "Tap1"))

n <- read_network(nd, de, hub, genes, neighbours = FALSE)

g <- graph_from_data_frame(n$edges, directed = TRUE, vertices = n$nodes)

# legend df for nodes
lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                     font.color = c("black", "white", "white", "white"),
                     shape = c( "ellipse"), 
                     color = c("orange", "darkred", "purple", "darkblue"),
                     title = "Informations", id = 1:4)

l <- layout_nicely(g)

pdf(file = paste0(basepath, "scripts/07_PlotsManuscript/plots_v2/12_Figure3_networkAbcTransporters.pdf"), 
    width = 10, height = 10)
plot(g, edge.arrow.size=.4, layout=l,
     vertex.size = 10,
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.9,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1.2,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
     edge.arrow.size = 0.8
)
legend(x=-1,y=-1,legend=c("DE", "hub", "no DE & no hub"),
       pch=21,
       col=c("orange", "darkred", "darkblue"), 
       pt.bg=c("orange", "darkred", "darkblue"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()



### CREATE NETWORK Tcf4 ####################################

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"

nd <- network_df("PFC")
de <- de_genes("PFC")
hub <- hub_genes("PFC")

genes <- reformat_genes(c("Tcf4"))

n <- read_network(nd, de, hub, genes)

g <- graph_from_data_frame(n$edges, directed = TRUE, vertices = n$nodes)

# legend df for nodes
lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                     font.color = c("black", "white", "white", "white"),
                     shape = c( "ellipse"), 
                     color = c("orange", "darkred", "purple", "darkblue"),
                     title = "Informations", id = 1:4)

l <- layout_nicely(g)

pdf(file = paste0(basepath, "scripts/07_PlotsManuscript/plots_v2/12_Figure5_networkTcf4.pdf"), 
    width = 10, height = 10)
plot(g, edge.arrow.size=.4, layout=l,
     vertex.size = 10,
     #vertex.label = ifelse(degree(g) > 4, V(g)$label, NA),
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.9,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1.2,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
     edge.arrow.size = 0.8
)
legend(x=-1,y=-1,legend=c("DE", "hub", "no DE & no hub"),
       pch=21,
       col=c("orange", "darkred", "darkblue"), 
       pt.bg=c("orange", "darkred", "darkblue"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


### CREATE NETWORK TRANSSYNAPTIC SIGNALING ------------

gene_str <- "CLSTN1:EPHB2:LRP8:NTNG1:RAP1A:SYT11:ATP1A2:ARF1:CACNB2:ITGB1:LRRC4C:SHANK2:GRIN2B:SYT1:SIPA1L1:PSEN1:VPS18:UNC13C:STX1B:CALB2:ABR:ARRB2:VAMP2:HAP1:CRHR1:MAPT:PPP1R9B:SYT4:PTPRS:CACNA1A:UNC13A:RAB3A:PRKCG:BRSK1:NRXN1:PLCB4:SNAP25:SRC:STAU1:KCNB1:APP:MAPK1:BCR:ADORA2A:CACNG2:SLC6A1:CSPG5:GNAI2:CACNA2D2:NLGN1:GRID2:PLK2:CAMK2A:GRIA1:CPLX2:RGS14:FLOT1:ITPR3:GRM4:HTR1B:GRIK2:GRM1:CAMK2B:STX1A:RELN:PRKAR2B:DGKI:ADRA1A:PTK2B:RIMS2:JAK2:STXBP1:CACNA1B:CASK:SYN1:SYP:NLGN3:FMR1"
genes <- stringr::str_split(gene_str, ":")[[1]]

nd <- network_df("vCA1")
de <- de_genes("vCA1")
hub <- hub_genes("vCA1")

genes <- reformat_genes(genes)

n <- read_network(nd, de, hub, genes, neighbours = FALSE)

g <- graph_from_data_frame(n$edges, directed = TRUE, vertices = n$nodes)

# legend df for nodes
lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                     font.color = c("black", "white", "white", "white"),
                     shape = c( "ellipse"), 
                     color = c("orange", "darkred", "#800080", "darkblue"),
                     title = "Informations", id = 1:4)

pdf(file = paste0(basepath, "scripts/07_PlotsManuscript/plots_v2/12_Figure4_networkTranssynapticSignaling_nicely.pdf"), 
    width = 10, height = 10)
plot(g, edge.arrow.size=.4, 
     # layout=layout_with_lgl,
     # layout=layout_with_dh,
     # layout=layout_with_gem,
     # layout=layout_with_graphopt,
     layout=layout_nicely,
     vertex.size = 8,
     #vertex.label = ifelse(degree(g) > 4, V(g)$label, NA),
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.9,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1.1,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
     edge.arrow.size = 0.8
)
legend(x=-1,y=-1,legend=c("DE", "hub", "DE & hub", "no DE & no hub"),
       pch=21,
       col=c("orange", "darkred", "#800080", "darkblue"), 
       pt.bg=c("orange", "darkred", "#800080", "darkblue"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

