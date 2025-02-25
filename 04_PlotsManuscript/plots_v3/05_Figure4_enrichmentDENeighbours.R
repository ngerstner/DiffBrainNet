##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GO enrichment plot for manuscript 
# vCA1 unique DE genes with neighbours enrichments


library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)

# 1. GO enrichment for vCA unique DE genes and their diff neighbours

data <- read_xlsx(path = "~/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure III_vCA1/network vCA1/Unique DEs and their diff neighbours_FUMA_gene2func67313/vCA1 Unique DE genes+diff neighbours_FUMA enrichments.xlsx",
                  sheet = "GO bp")

data <- arrange(data, desc(`[-LOG10 of FDR`))
data <- data[1:20,]

# remove the "GO_" from GO terms
data$GeneSet <- str_replace_all(data$GeneSet, '^GO_', ' ')
# remove the "_" from GO terms
data$GeneSet <- str_replace_all(data$GeneSet, "_", " ")
# only capitalize first letter per word
data$GeneSet <- str_to_title(data$GeneSet)
# of and to should start with lower case
data$GeneSet <- str_replace_all(data$GeneSet, "Of", "of")
data$GeneSet <- str_replace_all(data$GeneSet, "To", "to")
data$GeneSet <- str_replace_all(data$GeneSet, " In ", " in ")
data$GeneSet <- str_replace_all(data$GeneSet, " Into ", " into ")
# GeneSet as factor
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)


# color palette
colfunc <- colorRampPalette(c("#D45E60", "#0F5057"))

# dotplot enrichment
g <- ggplot(data) +
  geom_point(aes(y = reorder(GeneSet,`[-LOG10 of FDR`), 
                 x = `[-LOG10 of FDR`, 
                 color = `Genes ratio`, 
                 #color = `my adj p BH`, 
                 size = `OR[log2(a*d)-log2(b*c)]`)) +
  geom_vline(xintercept = -log10(0.05),linetype="dashed", color = "black")+
  scale_color_gradient(name = "% Genes Ratio",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "OddsRatio") +
  xlab("-log10 FDR") +
  ylab("GOterm") +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
ggsave(path= "~/ownCloud/DexStim_RNAseq_Mouse/scripts/07_PlotsManuscript/plots_v3", "05_Figure4_enrichmentDENeighbours_b.pdf", g, width = 8, height = 8)
