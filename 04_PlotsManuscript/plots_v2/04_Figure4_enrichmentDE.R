##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GO enrichment plot for manuscript 
# vCA1 unique DE genes enrichments

library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)

# 1. GO enrichment for vCA unique DE genes

data <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure III_vCA1/DE vCA1/FUMA_gene2func67310 (2)/vCA1 Unique DE genes_FUMA enrichments.xlsx",
                  sheet = "GO bp")

data <- arrange(data, desc(`[-log10 of FDR`))
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
# GeneSet as factor
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)


# color palette
colfunc <- colorRampPalette(c("#D45E60", "#0F5057"))

# dotplot enrichment
g <- ggplot(data) +
  geom_point(aes(y = reorder(GeneSet,`Genes ratio`), 
                 x = `Genes ratio`, 
                 color = `[-log10 of FDR`, 
                 size = `OR[log2(a*d)-log2(b*c)]`)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "OddsRatio") +
  xlab("% GeneRatio") +
  ylab("GOterm") +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
ggsave("04_Figure4_enrichmentDE.pdf", g, width = 8, height = 8)
