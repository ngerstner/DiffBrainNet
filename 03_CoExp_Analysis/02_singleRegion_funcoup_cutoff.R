##################################################
## Project: DexStim Mouse Brain
## Date: 03.12.2020
## Author: Nathalie
##################################################
# Decide on a beta cutoff for single region funcoup networks

library(data.table)
library(dplyr)
library(ggplot2)
library(igraph)
library(eulerr)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
region <- "PFC"
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

# 1. Read data
# 1a. Read co expression networks
dex0_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_singleRegion\\_",region,"\\_dex0\\_funcoup\\_parallel\\_.*\\.csv"),
                         full.names = TRUE)
dex1_files <- list.files(path = file.path(basepath, "tables/coExpression_kimono"), 
                         pattern = paste0("04\\_singleRegion\\_",region,"\\_dex1\\_funcoup\\_parallel\\_.*\\.csv"),
                         full.names = TRUE)

data_dex0 <- readFiles_concat(dex0_files)
data_dex1 <- readFiles_concat(dex1_files)


# 2. Join Base and Dex data frame
data <- inner_join(data_dex0, data_dex1, 
                   by = c("target", "predictor"),
                   suffix = c(".base", ".dex"))
dex_notBase <- anti_join(data_dex1, data_dex0, by = c("target", "predictor"))
head(data)

# Plot beta values
beta_data <- as.data.frame(rbind(cbind(rep("base", times = nrow(data)), data$beta_mean.base, data$relation.base),
                                 cbind(rep("dex", times = nrow(data)), data$beta_mean.dex, data$relation.dex)))
# ggplot(data = beta_data, aes(x = V1, y = V2, col = V3)) +
#   geom_boxplot() +
#   xlab("network") +
#   ylab("beta") +
#   theme(text = element_text(size=20))
# hist(data$beta_mean.base, breaks = 500) # histogram of beta values in base network
# png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
#                       region,"_funcoup_differential_betahistogramZoom.png"),
#     height = 400, width = 600)
# hist(data$beta_mean.base, breaks = 5000, xlim = c(-0.05,0.05),
#      xlab = "Beta", main = "", cex.lab=1.5, cex.axis = 1.2)
# dev.off()
# hist(data$beta_mean.dex, breaks = 500)  # histogram of beta values in dex network
# hist(data$beta_mean.dex, breaks = 5000, xlim = c(-0.05,0.05))


# 3. Remove interactions with very low r squared values & intercept & SVs
data <- data %>% 
  filter(overall_rsq.base >= rsquared_cutoff, overall_rsq.dex >= rsquared_cutoff) %>%
  filter(predictor != '(Intercept)')
data <- data[!startsWith(data$predictor, "SV"),]

# 4. Calculate z scores for interactions that are left
data <- mutate(data, z = z_score(beta_mean.dex,beta_mean.base, beta_stderr.dex, beta_stderr.base))
hist(data$beta_mean.base)

# 5. Cutoff optimization
poss_cutoff <- c(0.000001, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.1)
nodebetweenness_list <- list()
for(beta_cutoff in poss_cutoff) {
  
  print(beta_cutoff)
  
  # 5a. Keep only interactions that have at least in one network a beta value > 0.05
  data <- data %>%
    mutate(diff = (abs(beta_mean.base) > beta_cutoff | abs(beta_mean.dex) > beta_cutoff))
  data_diff1 <- data %>%
    filter(diff)
  data_diff1$p_diff <- 2*pnorm(-abs(data_diff1$z))
  data_diff1$p_adj <- p.adjust(data_diff1$p_diff, method = "fdr")
  
  # 5b. Create network corresponding to beta cutoff
  data_diff <- data_diff1 %>%
    filter(p_adj <= padj_cutoff)
  head(data_diff[,c("beta_mean.base", "beta_stderr.base", "beta_mean.dex", "beta_stderr.dex", "z", "p_adj")], 20)
  
  # Nodes in network
  node_vec <- unique(c(data_diff$target, data_diff$predictor))
  
  # Find modules
  relations <- data.frame(from=data_diff$target,
                          to=data_diff$predictor,
                          value=data_diff$z,
                          performance=data_diff$p_adj) 
  g_diff <- graph_from_data_frame(relations, directed=FALSE, vertices=node_vec)
  
  # Calculate nodedegree
  nodedegree <- igraph::degree(g_diff)
  nodedegree <- sort(nodedegree, decreasing = TRUE)
  
  # Calculate nodebetweenness
  nodebetweenness <- betweenness(g_diff, directed = FALSE)  # node betweenness: number of shortest paths going through a node
  nodebetweenness <- sort(nodebetweenness, decreasing = TRUE)
  nodebetweenness_list[[as.character(beta_cutoff)]] <- nodebetweenness
  # hist(nodebetweenness)
}


# Compare top genes between different cutoffs
for (i in 1:(length(nodebetweenness_list)-1)){
  
  intersect_top <- intersect(names(nodebetweenness_list[[i]][1:100]), names(nodebetweenness_list[[i+1]][1:100]))
  print(length(intersect_top))
  match_pos <- match(names(nodebetweenness_list[[i]])[names(nodebetweenness_list[[i]]) %in% intersect_top],
        names(nodebetweenness_list[[i+1]])[names(nodebetweenness_list[[i+1]]) %in% intersect_top])
  print(match_pos)
  
  # Rank plot to compare ranks of top genes
  # Comparison between i and i+1
  png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                            region,"_funcoup_cutoffSelection_rankPlot_betacutoffs",poss_cutoff[i],"-",poss_cutoff[i+1],".png"),
      width = 700, height = 700)
  plotRanks(names(nodebetweenness_list[[i]][1:100]), names(nodebetweenness_list[[i+1]][1:100]),
            labels = FALSE)
  dev.off()
  # Comparison between 1 and i
  if (i > 1){
    png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                          region,"_funcoup_cutoffSelection_rankPlot_betacutoffs",poss_cutoff[1],"-",poss_cutoff[i],".png"),
        width = 700, height = 700)
    plotRanks(names(nodebetweenness_list[[1]][1:100]), names(nodebetweenness_list[[i]][1:100]),
              labels = FALSE)
    dev.off()
  }
  
  # Venn diagram (euler plot)
  # Comparison between i and i+1
  list1 <- list(names(nodebetweenness_list[[i]][1:100]), 
                names(nodebetweenness_list[[i+1]][1:100]))
  names(list1) <- c(paste("Beta cutoff", poss_cutoff[i]),
                    paste("Beta cutoff", poss_cutoff[i+1]))
  png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                        region,"_funcoup_cutoffSelection_vennEuler_betacutoffs",poss_cutoff[i],"-",poss_cutoff[i+1],".png"),
      width = 700, height = 700)
  print(plot(euler(list1, shape = "ellipse"), 
             labels = list(cex = 1.5), quantities = list(cex = 1.5)))
  dev.off()
  
  # Comparison between i and 1
  if(i > 1){
    list1 <- list(names(nodebetweenness_list[[1]][1:100]), 
                  names(nodebetweenness_list[[i]][1:100]))
    names(list1) <- c(paste("Beta cutoff", poss_cutoff[1]),
                      paste("Beta cutoff", poss_cutoff[i]))
    png(filename = paste0(basepath, "figures/02_CoExp_Kimono/03_AnalysisFuncoup/singleRegions/",
                          region,"_funcoup_cutoffSelection_vennEuler_betacutoffs",poss_cutoff[1],"-",poss_cutoff[i],".png"),
        width = 700, height = 700)
    print(plot(euler(list1, shape = "ellipse"), 
               labels = list(cex = 1.5), quantities = list(cex = 1.5)))
  dev.off()
  }
}

# I WOULD DECIDE ON 0.001 from plots, but top genes have too many connections with 0.001
# --> probably 0.01