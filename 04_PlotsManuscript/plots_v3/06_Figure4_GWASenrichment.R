##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GWAS enrichment plot for manuscript 
# vCA1

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

# 1. Pathway enrichment for Abcd1 and neighbours

data <- read_xlsx(path = "~/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure III_vCA1/network vCA1/Unique DEs and their diff neighbours_FUMA_gene2func67313/vCA1 Unique DE genes+diff neighbours_FUMA enrichments.xlsx",
                  sheet = "GWAS")

data <- arrange(data, desc(`[-log10 of FDR`))
data <- data[data$`my adj p BH` <= 0.05,]
data <- data[!is.na(data$`my adj p BH`),]

# GeneSet as factor 
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)

# color palette
colfunc <- colorRampPalette(c("#D45E60", "#0F5057"))

# dotplot enrichment
g <- ggplot(data) +
  geom_point(aes(y = reorder(GeneSet,`[-log10 of FDR`), 
                 x = `[-log10 of FDR`, 
                 color = `Genes ratio`, 
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
ggsave(path= "~/ownCloud/DexStim_RNAseq_Mouse/scripts/07_PlotsManuscript/plots_v3", "06_Figure4_GWASenrichment.pdf", g, width = 8, height = 5)
