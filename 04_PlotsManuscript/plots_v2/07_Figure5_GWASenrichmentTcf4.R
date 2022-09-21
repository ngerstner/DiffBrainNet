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

data <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure IV_GWAS/TCF4/FUMA_diffTcf4PFCnetwork_gene2func68788/FUMA_diffPFCTcf4.xlsx",
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
  geom_point(aes(y = reorder(GeneSet,`Genes ratio`), 
                 x = `Genes ratio`, 
                 color = `[-log10]`, 
                 size = `OR[log2(a*d)-log2(b*c)]`)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "OddsRatio") +
  xlab("% GeneRatio") +
  ylab("GOterm") +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
ggsave("07_Figure5_GWASenrichmentTcf4.pdf", g, width = 8, height = 6)
