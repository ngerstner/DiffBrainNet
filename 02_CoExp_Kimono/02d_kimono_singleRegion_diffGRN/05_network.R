##################################################
## Project: DexStim Mouse Brain
## Date: 29.09.2020
## Author: Nathalie
##################################################
# Analyze kimono output and create network
# Separate on dex and baseline

library(data.table)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)

reg <- "AMY"
d <- 0

basepath <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
kimono_results <- paste0(basepath,"tables/coExpression_kimono/04_singleRegion_",reg,"_dex",d,"_funcoup.csv")
GR_genes <- paste0(basepath,"data/kimono_input/63genes_ZimmermannPaper.csv")

# 1. Load Kimono results and filter them
d <- fread(kimono_results)
#dtmoni <- d[value!=0,]; nrow(dtmoni) 
dtmoni <- d
network <- dtmoni %>% #filter((value > 0.001) | (value < (-0.001))) %>%
  #filter(performance > 0.01) %>%
  filter(predictor != '(Intercept)') %>% setDT

png(filename = paste0(basepath, "figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_allbetas.png"))
hist(abs(network$value), breaks = 500, 
     main = paste("Distribution of all beta values in", reg),
     xlab = "Beta",xlim = c(0,0.1))
abline(v = mean(abs(network$value)), lty = 3)
dev.off()

hist(table(network$target))


# 2. Inspect associations
dtmoni$target %>%  unique %>%  length
dtmoni[relation=="mrna_bio", predictor] %>% unique

ggplot(network, aes(relation)) + geom_bar(fill="lightblue")

network[,.N, by=relation]
network[, .N, by=predictor][order(-N)][1:100]

nrow(network)/nrow(dtmoni) # retained fraction after filtering

paste("Number of unique genes:", 
      length(unique(network$target)))
paste("Number of unique predictors:", 
      length(unique(network$predictor)))

# unique predictors
unique(network, by=c("relation", "predictor")) %>% .[,.N, by=relation]


# 3. Create network
library(purrr)
library(igraph)
unique(network$relation)

actors<-unique(c(network$predictor,network$target))

relations <- data.frame(from=network$predictor,
                        to=network$target,
                        value=network$value,
                        performance=network$performance) 

# network
g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
print(g)


# 4. Plot network
plot_network <- function(subnet){
  subnet <- data.frame(lapply(subnet, as.character), stringsAsFactors=FALSE)
  subnet$target <- gsub("_", ".", subnet$target)
  subnet$predictor <- gsub("_", ".", subnet$predictor)
  subnet$relation <- gsub("_",".", subnet$relation)
  
  lable <- unique(c(paste0("mrna_",subnet[,1]),paste0(subnet$relation,"_",subnet[,2])))
  
  actors <- data.frame(name=do.call(rbind, strsplit(lable,"_") )[,2],
                       omic= do.call(rbind, strsplit(lable,"_") )[,1]
  )
  
  
  relations <- data.frame(from=subnet$predictor,
                          to=subnet$target,
                          value=subnet$value,
                          performance=subnet$performance) 
  
  actors <- unique(actors) %>%  setDT
  actors[omic=="prior.mrna", omic:="mrna"]
  actors <- unique(actors)
  g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  deg2 <- igraph::degree(g, mode="all")
  
  # plot(g, vertex.label.color="black", vertex.size=10, vertex.label=NA , vertex.label.dist=1.5)
  # library(qgraph)
  
  e <- get.edgelist(g,names=FALSE)
  # l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),  area=1000*(vcount(g)^2),repulse.rad=100+(vcount(g)^30))
  
  library(RColorBrewer)
  
  actors[, edges:=igraph::degree(g)]
  actors[, mynames:=name]
  
  
  actors$id <- 1:nrow(actors)
  actors <-actors[order(id)]
  summary(actors$edges)
  
  
  # edge_threshold <- (as.numeric(sub('.*:', '', summary(actors$edges)[5])) +as.numeric(sub('.*:', '', summary(actors$edges)[6])))/2
  
  #edge_threshold<-500
  #actors[edges>edge_threshold,]
  #actors[,.N, by=edges,][order(edges)]
  #actors[edges < edge_threshold  , mynames:=NA]
  
  colrs <-c(brewer.pal(4, "Set2")[c(1:4)])
  mycol=colrs[2:3]
  actors[omic=="mrna", mycolor:=colrs[2]]
  actors[omic=="mrna.bio", mycolor:=colrs[3]]
  plot(g,
       layout=layout.fruchterman.reingold(g),
       vertex.frame.color= adjustcolor("black", .4)	, 
       vertex.size=1+(log(deg2)*2),
       vertex.color=actors$mycolor ,
       edge.color =  adjustcolor("grey", .8),
       edge.curved=.1,
       vertex.label = actors$mynames,
       vertex.label.color= adjustcolor("black", .8),
       vertex.label.cex = 0.7,    asp = 1 ,
       # vertex.label.family = "Times",
       edge.width=E(g)$weight*400,
       main = paste("Number of edges:", nrow(subnet)))
  legend(x=-1, y=-0.5,c("gene", "covariable"),
         pch=21, col="#777777", pt.bg=mycol, pt.cex=2, cex=.8, bty="n", ncol=1)
}

# Function to change all ensembl ids in network to gene symbols
ensemblToSymbol <- function(net){
  net$target <- mapIds(org.Mm.eg.db, keys = net$target, 
                       keytype = "ENSEMBL", column="SYMBOL")
  net$predictor[grepl("ENSM", net$predictor)] <- mapIds(org.Mm.eg.db, keys = net$predictor[grepl("ENSM", net$predictor)], 
                                                        keytype = "ENSEMBL", column="SYMBOL")
  return(net)
}

# whole network
#plot_network(network)

# dex network
dex <-network[predictor=="Dex",]
#plot_network(dex)

# Read 63 genes from Zimmermann/Arloth Paper
GRgenes <- fread(GR_genes)
# Plot the 63 genes and all their connections
net_GR <- network[predictor %in% GRgenes$Ensembl | target %in% GRgenes$Ensembl,]
#png(filename = paste0(basepath,"figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_63withpredictors.png"))
plot_network(net_GR)
#dev.off()

# Plot the 63 genes, regions and dex, and the connections between 
net_GR <- network[target %in% GRgenes$Ensembl,]
# net_GR <- net_GR[predictor %in% GRgenes$Ensembl, ]
net_GR <- net_GR[predictor %in% GRgenes$Ensembl , ]
plot_network(net_GR)

png(filename = paste0(basepath, "figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_63betas.png"))
hist(abs(net_GR$value), breaks = 20, 
     main = paste("Distribution of beta values in", reg),
     xlab = "Beta")
abline(v = mean(abs(net_GR$value)), lty = 3)
dev.off()

png(filename = paste0(basepath, "figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_63performance.png"))
hist(net_GR$performance, 
     main = paste("Distribution of performance values in", reg),
     xlab = "Performance")
dev.off()

net_GR <- ensemblToSymbol(net_GR)
png(filename = paste0(basepath,"figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_63wofilter.png"),
    width = 900, height = 900)
plot_network(net_GR)
dev.off()

net_GR <- net_GR[abs(value) >= 0.001,]
png(filename = paste0(basepath,"figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_63filter0.001.png"),
    width = 900, height = 900)
plot_network(net_GR)
dev.off()

# Plot with node degress
degrees <- count(network, target) %>%
  mutate(network = "complete")
boxplot(degrees$n)
degrees_prior <- count(network[network$relation == "prior_mrna",], target) %>%
  mutate(network = "prior")
boxplot(degrees_prior$n)
degrees_bio <- count(network[network$relation == "mrna_bio",], target) %>% 
  mutate(network = "biological")
boxplot(degrees_bio$n)

degrees_comb <- rbind(degrees, degrees_prior, degrees_bio)
degrees_comb$network <- factor(degrees_comb$network, levels = c("complete", "prior", "biological"))

ggplot(degrees_comb, aes(x = network, y = n, fill = network)) +
  geom_boxplot() +
  ylab("node degree")
ggsave(filename = paste0(basepath,"figures/02_CoExp_Kimono/05_singleRegion_",reg,"_funcoup_nodeDegrees_wofilter.png"))



# Compare number of edges between 63 genes with number of edges between randomly selected 63 genes
mrna <- d[d$relation == "prior_mrna",]
genes <- unique(d$predictor)
