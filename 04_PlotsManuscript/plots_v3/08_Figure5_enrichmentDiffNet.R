##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GO enrichment plot for manuscript 
# Tcf4 differential network in PFC

library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)

# 1. GO enrichment for vCA unique DE genes and their diff neighbours

data <- read_xlsx(path = "~/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure IV_GWAS/TCF4/FUMA_diffTcf4PFCnetwork_gene2func68788/FUMA_diffPFCTcf4.xlsx",
                  sheet = "GO bp")

data <- arrange(data, desc(`[-log10]`))
data <- data[1:25,]

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
# make some terms shorter with abbreviations
data$GeneSet <- str_replace_all(data$GeneSet, "Positive Regulation", "Pos. Reg.")
data$GeneSet <- str_replace_all(data$GeneSet, "Negative Regulation", "Neg. Reg.")
# Ii from "Polymerase II" back to upper case
data$GeneSet <- str_replace_all(data$GeneSet, "Ii", "II")
# GeneSet as factor
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)


# color palette
colfunc <- colorRampPalette(c("#D45E60", "#0F5057"))

# dotplot enrichment
g <- ggplot(data) +
  geom_point(aes(y = reorder(GeneSet,`[-log10]`), 
                 x = `[-log10]`, 
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
ggsave(path= "~/ownCloud/DexStim_RNAseq_Mouse/scripts/07_PlotsManuscript/plots_v3", "08_Figure5_enrichmentDiffNet_NG.pdf", g, width = 8, height = 8)
