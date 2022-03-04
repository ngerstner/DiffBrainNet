##################################################
## Project: DexStim Mouse Brain
## Date: 29.09.2020
## Author: Nathalie
##################################################
# Parse FunCoup mouse data as input for kimono

library(data.table)
library(dplyr)
library(org.Mm.eg.db)
library(igraph)

basepath <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/"
funcoup_file <- file.path(basepath, "data/kimono_input/FC5.0_M.musculus_compact.gz")
output_file <- file.path(basepath,"data/kimono_input/prior_expr_funcoup_mm.csv")

# 1. Read file
funcoup_mus <- fread(funcoup_file, col.names = c("PFC", "FBS", "Gene_A", "Gene_B")) #%>%
  #filter(PFC >= 0.4)


# 2. Write interactions with ENSEMBL IDs to file
funcoup_ens <- funcoup_mus %>%
  dplyr::select('Gene_A', 'Gene_B') 
fwrite(funcoup_ens, file = output_file)


# 3. Plot some statistics
funcoup_app <- c(funcoup_ens$Gene_A, funcoup_ens$Gene_B)
funcoup_unique <- unique(funcoup_app)
hist(table(funcoup_app))
max(table(funcoup_app))


# 4. Read 63 genes from Zimmermann/Arloth Paper
# And check if they interact in GeneMANIA
GR_genes <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/63genes_ZimmermannPaper.csv"
GRgenes <- fread(GR_genes)
funcoup_GR <- funcoup_mus %>% filter(Gene_A %in% GRgenes$Ensembl & Gene_B %in% GRgenes$Ensembl)


# Function to change all ensembl ids in network to gene symbols
ensemblToSymbol <- function(net){
  net$Gene_A <- mapIds(org.Mm.eg.db, keys = net$Gene_A, 
                       keytype = "ENSEMBL", column="SYMBOL")
  net$Gene_B <- mapIds(org.Mm.eg.db, keys = net$Gene_B,
                       keytype = "ENSEMBL", column="SYMBOL")
  return(net)
}

funcoup_GR <- ensemblToSymbol(funcoup_GR)
funcoup_GR <- funcoup_GR[,c("Gene_A", "Gene_B")]
df.g <- graph.data.frame(d = funcoup_GR, directed = FALSE)

png(filename = paste0(basepath,"figures/02_CoExp_Kimono/05_funcoup_63genes.png"), width = 900, height = 900)
plot(df.g, vertex.label = V(df.g)$name,
     layout=layout.fruchterman.reingold(df.g),
     vertex.label.color= adjustcolor("black", .8),
     vertex.label.cex = 0.7,
     vertex.size=10)
dev.off()
