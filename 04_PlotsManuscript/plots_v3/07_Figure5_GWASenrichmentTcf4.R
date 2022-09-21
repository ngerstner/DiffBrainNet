##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GWAS enrichment plot for manuscript 
# Tcf4 neighbours?

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

# 1. Pathway enrichment for Abcd1 and neighbours

data <- read_xlsx(path = "~/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure IV_GWAS/TCF4/FUMA_diffTcf4PFCnetwork_gene2func68788/FUMA_diffPFCTcf4.xlsx",
                  sheet = "GWAS")

data <- arrange(data, desc(`[-log10]`))
data$`my adj p` <- as.numeric(data$`my adj p`)
data <- data[data$`my adj p` <= 0.05,]
data <- data[!is.na(data$`my adj p`),]

# GeneSet as factor 
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)


# color palette
colfunc <- colorRampPalette(c("#D45E60", "#0F5057"))

# dotplot enrichment
g <- ggplot(data) +
  geom_point(aes(y = reorder(GeneSet,`[-log10]`), 
                 x = `[-log10]`, 
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
ggsave(path= "~/ownCloud/DexStim_RNAseq_Mouse/scripts/07_PlotsManuscript/plots_v3", "07_Figure5_GWASenrichmentTcf4.pdf", g, width = 8, height = 6)
