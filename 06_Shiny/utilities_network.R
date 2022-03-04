### FUNCTIONS -------------------------------------


# remove duplicated edges from network dataframe
simplify_network_df <- function(network_dt){
  
  network_dt <- network_dt %>%
    dplyr::filter(from != to) #%>%
    # dplyr::mutate(from_sort = pmin(from, to), to_sort = pmax(from, to)) %>%
    # dplyr::distinct(from_sort, to_sort, .keep_all = TRUE) %>%
    # dplyr::select(!c(from_sort, to_sort))
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


# function to generate baseline/diff network, single regions
read_network <- function(data, de_genes, hub_genes, gene, mode, id_type, neighbours){
  
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
  print(nrow(relations))
  
  # Nodes
  # Generate dataframe for node layout
  node_vec <- unique(c(relations$from, relations$to))
  if (id_type == "Gene Symbol"){
    labels <- mapIds(org.Mm.eg.db, keys = node_vec,
                     column = "SYMBOL", keytype = "ENSEMBL")
  } else {
    labels <- node_vec
  }
  
  # set colours according to DE status of gene
  col <- rep("darkblue", length(node_vec))
  col[node_vec %in% de_genes] <- "orange"
  col[node_vec %in% hub_genes] <- "darkred"
  col[node_vec %in% de_genes & node_vec %in% hub_genes] <- "purple"
  nodes <- data.table(id = node_vec,
                      label = labels,
                      color = col)
  
  # nodes data.frame for legend
  lnodes <- data.frame(label = c("DE", "hub", "DE & hub", "no DE & no hub"),
                       font.color = c("black", "white", "white", "white"),
                       shape = c( "ellipse"), 
                       color = c("orange", "darkred", "purple", "darkblue"),
                       title = "Informations", id = 1:4)
                      
  return(list("edges" = relations, "nodes" = nodes, "ledges" = ledges, "lnodes" = lnodes))
  
}


# function to generate baseline/diff network, multiple regions
read_network_multi <- function(regions, gene_unsplit, mode, id_type){
  
  list_reg <- list()
  
  for (reg in regions){
    
    # Read data
    data <-
      fread(
        file = paste0(
          basepath,
          "tables/coExpression_kimono/03_AnalysisFuncoup/",
          "04_singleRegion_",
          reg,
          "_filtered_diffNetwork.csv"
        )
      ) %>%
      dplyr::select(target, predictor, beta_mean.base, beta_mean.dex, z, p_adj) 
    
    cols <- colnames(data)[3:6]
    data <- data %>%
      dplyr::rename_with(.fn = ~paste0(., ".", reg), .cols = all_of(cols) )
    
    list_reg[[reg]] <- data
    
  }
  
  data_joined <- list_reg %>% 
    purrr::reduce(full_join, by = c("target","predictor"))
  
  # # Read data
  # data <- fread(file = paste0(basepath, "tables/coExpression_kimono/03_AnalysisFuncoup/",
  #                             "04_singleRegion_",region,"_filtered_", mode, "Network.csv"))
  gene <- stringr::str_split(gene_unsplit, ",\\s?")[[1]]
  print(gene)
  if (!startsWith(gene_unsplit, "ENSMUSG")){
    gene <- mapIds(org.Mm.eg.db, keys = gene, column = "ENSEMBL", keytype = "SYMBOL")
  }
  
  # Edges
  c <- rep("#db9512", nrow(data))
  c[data$z > 0] <- "#128bdb"
  relations <- data.table(from=data$target,
                          to=data$predictor,
                          z = data$z,
                          p_adj = data$p_adj,
                          color = c)
  # print(relations)
  relations <- simplify_network_df(relations)
  relations <- subset_network(relations, gene)
  print(nrow(relations))
  
  # Nodes
  node_vec <- unique(c(relations$from, relations$to))
  if (id_type == "Gene Symbol"){
    labels <- mapIds(org.Mm.eg.db, keys = node_vec,
                     column = "SYMBOL", keytype = "ENSEMBL")
  } else {
    labels <- node_vec
  }
  nodes <- data.table(id = node_vec,
                      label = labels)
  # title = node_vec)
  return(list("edges" = relations, "nodes" = nodes))
  
}


# function to generate network stats
network_stats <- function(df){
  
  # Initialize data frame for stats
  stats <- data.frame(#network = character(),
                      type = character(),
                      value = numeric())
  
  stats <- rbind(stats, list("type" = "edges", "value" = nr_edges_base))
}
