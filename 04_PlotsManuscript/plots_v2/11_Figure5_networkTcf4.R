##################################################
## Project: DexStim Mouse Brain
## Date: 03.03.2022
## Author: Nathalie
##################################################
# Tcf4 network in vDG and vCA1 with igraph to pdf


library(igraph)
library(data.table)
library(dplyr)
library(org.Mm.eg.db)
library(RColorBrewer)


### FUNCTIONS ####################################

network_data <- function(regions, network_type = "diff"){
  
  # Read data tables from all required regions
  list_network <- list()
  
  for (reg in regions){
    file_path <- paste0(
      basepath,
      "tables/coExpression_kimono/03_AnalysisFuncoup/",
      "04_singleRegion_",
      reg,
      "_filtered_",
      network_type,
      "Network.csv"
    )
    # Read data
    if (network_type == "diff") {
      data <-
        fread(file = file_path) %>%
        dplyr::select(target, predictor, beta_mean.base, beta_mean.dex, z, p_adj)
      #dplyr::select(target, predictor, z) %>%
      #dplyr::rename_with(.fn = ~ reg, .cols = z)
      
      cols <- colnames(data)[3:6]
      data <- data %>%
        dplyr::rename_with(.fn = ~paste0(., ".", reg), .cols = all_of(cols) )
    } else {
      data <-
        fread(file = file_path) %>%
        dplyr::select(target, predictor, beta_mean)
      # cols <- colnames(data)[3]
      data <- data %>%
        dplyr::rename_with(.fn = ~ reg, .cols = beta_mean)
      # dplyr::rename_with(.fn = ~paste0(., ".", reg), .cols = all_of(cols) )
    }
    list_network[[reg]] <- data
  }
  
  data_joined <- list_network %>% 
    purrr::reduce(full_join, by = c("target","predictor"))
  
  return(data_joined)
}


# DE genes of required brain regions
de_genes <- function(regions){
  
  list_de <- list()
  # Read and subset DE genes for coloring
  for (reg in regions){
    df_de <- fread(paste0(basepath, "tables/02_",
                          reg,"_deseq2_Dex_1_vs_0_lfcShrink.txt"))
    na_indices <- which(is.na(df_de$padj))
    df_de$padj[na_indices] <- 1
    df_de <- df_de[df_de$padj <= 0.1,]
    
    df_de <- df_de %>%
      dplyr::select(Ensembl_ID, padj) %>%
      dplyr::rename_with(.fn = ~ reg, .cols = padj)
    
    list_de[[reg]] <- df_de
  }
  
  de_joined <- list_de %>% 
    purrr::reduce(full_join, by = c("Ensembl_ID"))
  
  return(de_joined)
}


# hub genes of required brain regions
hub_genes <- function(regions){
  
  list_hub <- list()
  # Read and subset hub genes for coloring
  for (reg in regions) {
    df_hub <-
      fread(
        paste0(
          basepath,
          "tables/coExpression_kimono/03_AnalysisFuncoup/04_",
          reg,
          "_funcoup_differential",
          "_nodebetweennessNorm_betacutoff0.01.csv"
        )
      ) %>%
      filter(nodebetweenness_norm >= 1) %>%
      dplyr::select(ensembl_id, nodebetweenness_norm) %>%
      dplyr::rename_with(.fn = ~ reg, .cols = nodebetweenness_norm)
    
    list_hub[[reg]] <- df_hub
  }
  
  hub_joined <- list_hub %>% 
    purrr::reduce(full_join, by = c("ensembl_id"))
  
  return(hub_joined)
}


# bring input genes to correct format
reformat_genes <- function(list_genes){
  if (!startsWith(list_genes[1], "ENSMUSG")){
    format_gene <- sapply(list_genes, stringr::str_to_title)
    format_gene <- mapIds(org.Mm.eg.db, keys = format_gene, column = "ENSEMBL", keytype = "SYMBOL")
    ids_na <- names(format_gene)[which(is.na(format_gene))]
    #print(ids_na)
    if (length(ids_na) > 0) {
      showNotification(paste("No Ensembl ID found for following genes: ",
                             ids_na), type = "message", duration = 5)
    }
    format_gene <- format_gene[which(!is.na(format_gene))]
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


# create network function
network <- function(data, de_genes, hub_genes, input_genes, mode = "diff", 
                    id_type = "Gene Symbol", neighbours = TRUE){
  
  # input genes
  print(input_genes)
  
  # filter data
  data <- data %>%
    dplyr::rename(from = predictor, to = target) %>%
    dplyr::filter(from != to) 
  data <- subset_network(data, input_genes, neighbours)
  
  # color palettes for edges and nodes
  # edge palette --> assumes that data stores only target, predictor and z scores
  # edge_colors <- brewer.pal(ncol(data)-1, "Blues")[2:(ncol(data)-1)]
  
  # Edges
  if(mode == "baseline" | mode == "treatment"){
    # color palettes for edges and nodes
    edge_colors <- brewer.pal(ncol(data)-1, "Blues")[2:(ncol(data)-1)]
    
    # count regions per diff edge that are not na
    data$edge_reg <- apply(data[,3:ncol(data)], 1, function(x) names(which(!is.na(x))) )
    data$count_reg <- sapply(data$edge_reg, length)
    data$title <- paste0("<p><b>Connection in: </b>", sapply(data$edge_reg, paste, collapse = ", "),
                         "</p>")
    # assign color to edges according to number of regions
    data$c <- sapply(data$count_reg, function(x) edge_colors[x])
    
    print(head(data))
    # df for edges
    relations <- data.table(from=data$from,
                            to=data$to,
                            #beta = data$beta_mean,
                            color = c,
                            title = data$title)
    # df for edge legend
    ledges <- data.frame(color = edge_colors,
                         label = 1:(ncol(data)-6), 
                         font.align = "top")
  } else {
    # remove columns
    data <- data %>% dplyr::select(-contains("beta"), -contains("p_adj"))
    
    # color palettes for edges and nodes
    # edge palette --> assumes that data stores only target, predictor and z scores
    edge_colors <- brewer.pal(ncol(data)-1, "Blues")[2:(ncol(data)-1)]
    
    # count regions per diff edge that are not na
    data$edge_reg <- apply(data[,3:ncol(data)], 1, function(x) names(which(!is.na(x))) )
    data$count_reg <- sapply(data$edge_reg, length)
    data$title <- paste0("<p><b>Diff. co-expressed in: </b>", sapply(data$edge_reg, paste, collapse = ", "),
                         "</p>")
    # assign color to edges according to number of regions
    data$c <- sapply(data$count_reg, function(x) edge_colors[x])
    
    print(head(data))
    # df for edges
    relations <- data.table(from=data$from,
                            to=data$to,
                            #z = data$z,
                            #p_adj = data$p_adj,
                            color = data$c,
                            title = data$title)
    # df for edge legend
    ledges <- data.frame(color = edge_colors,
                         label = 1:(ncol(data)-6), 
                         font.align = "top")
  }
  print(nrow(relations))
  
  
  # Nodes
  # get unique nodes with correct id
  nodes <- data.frame("id" = unique(union(
    c(relations$from, relations$to),
    input_genes
  )), stringsAsFactors = FALSE)
  if (id_type == "Gene Symbol"){
    print(nodes$id)
    nodes$label <- mapIds(org.Mm.eg.db, keys = nodes$id,
                          column = "SYMBOL", keytype = "ENSEMBL")
  } else {
    nodes$label <- nodes$id
  }
  
  # count regions where gene is DE or/and hub
  nodes_de <- left_join(x = nodes, y = de_genes,
                        by = c("id" = "Ensembl_ID"))
  nodes_hub <- left_join(x = nodes, y = hub_genes,
                         by = c("id" = "ensembl_id"))
  
  nodes$de_reg <- apply(nodes_de[,3:ncol(nodes_de)], 1, 
                        function(x) names(which(!is.na(x))) )
  nodes$de_count <- sapply(nodes$de_reg, length)
  
  nodes$hub_reg <- apply(nodes_hub[,3:ncol(nodes_hub)], 1, 
                         function(x) names(which(!is.na(x))) )
  nodes$hub_count <- sapply(nodes$hub_reg, length)
  
  # set color according to DE/hub status
  nodes$color <- rep("darkblue", nrow(nodes_de))
  nodes$color[nodes$de_count > 0] <- "orange"
  nodes$color[nodes$hub_count > 0] <- "darkred"
  nodes$color[nodes$de_count > 0 & nodes$hub_count > 0] <- "purple"
  
  #nodes$opacity <- nodes$de_count + nodes$hub_count
  #nodes$opacity <- nodes$opacity/max(nodes$opacity)
  
  nodes$title <- paste0("<p><i>",nodes$label,"</i>",
                        "<br><b>DE in regions: </b>", sapply(nodes$de_reg, paste, collapse = ", "),
                        "<br><b>Hub in regions: </b>", sapply(nodes$hub_reg, paste, collapse = ", "),"</p>")
  
  print(head(nodes))
  
  # df for node data frame
  lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                       shape = c( "ellipse"), color = c("orange", "darkred", "purple", "darkblue"),
                       title = "Informations", id = 1:4)
  
  
  # return(list("edges" = relations, "nodes" = nodes, "ledges" = ledges, "lnodes" = lnodes))
  return(list("edges" = relations, "nodes" = nodes, "ledges" = ledges))
  
}



### CREATE NETWORK #######################################

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"

nd <- network_data(c("dDG", "vDG"))
de <- de_genes(c("dDG", "vDG"))
hub <- hub_genes(c("dDG", "vDG"))

genes <- reformat_genes(c("Tcf4"))

n <- network(nd, de, hub, genes)

g <- graph_from_data_frame(n$edges, directed = TRUE, vertices = n$nodes)

# legend df for nodes
lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                     font.color = c("black", "white", "white", "white"),
                     shape = c( "ellipse"), 
                     color = c("orange", "darkred", "purple", "darkblue"),
                     title = "Informations", id = 1:4)

l <- layout_nicely(g)

pdf(file = "11_Figure5_networkTcf4_partlyLabels.pdf", width = 10, height = 10)
plot(g, edge.arrow.size=.4, layout=l,
     vertex.size = 6,
     vertex.label = ifelse(degree(g) > 4, V(g)$label, NA),
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.7,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
)
legend(x=-1,y=-1,legend=c("DE", "hub", "no DE & no hub"),
       pch=21,
       col=c("orange", "darkred", "darkblue"), 
       pt.bg=c("orange", "darkred", "darkblue"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

pdf(file = "11_Figure5_networkTcf4_woLabels.pdf", width = 10, height = 10)
plot(g, edge.arrow.size=.4, layout=l,
     vertex.size = 6,
     vertex.label = ifelse(degree(g) > 50, V(g)$label, NA),
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.7,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
)
dev.off()

pdf(file = "11_Figure5_networkTcf4_wLabels.pdf", width = 10, height = 10)
plot(g, edge.arrow.size=.4, layout=l,
     vertex.size = 6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
     #vertex.label.font=c(1,2,3,4),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.7,                 # Font size (multiplication factor, device-dependent)
     vertex.label.dist=1,                           # Distance between the label and the vertex
     vertex.label.degree=80,                        # The position of the label in relation to the vertex (use pi))
)
dev.off()
