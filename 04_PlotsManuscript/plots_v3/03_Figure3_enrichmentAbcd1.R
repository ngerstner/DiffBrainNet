##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# Pathway enrichment dotplot for manuscript 
# Abcd1 and neighbourhood in diff network of PFC

library(readxl)

# 1. Pathway enrichment for Abcd1 and neighbours

data <- read_xlsx(path = "~/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/Abcd1/Abcd1 gene neighbourhood_FUMA_gene2func67210/Abcd1_diffnetwork_FUMA.xlsx",
                  sheet = "KEGG AND REACTOME")

data <- arrange(data, desc(`[-log10FDR]`))
data <- data[1:20,]


# remove the "_" from GO terms
data$GeneSet <- str_replace_all(data$GeneSet, "_", " ")
# only capitalize first letter per word
data$GeneSet <- str_to_title(data$GeneSet)
# of and to should start with lower case
data$GeneSet <- str_replace_all(data$GeneSet, "Of", "of")
data$GeneSet <- str_replace_all(data$GeneSet, "To", "to")
data$GeneSet <- str_replace_all(data$GeneSet, " In ", " in ")
data$GeneSet <- str_replace_all(data$GeneSet, " Into ", " into ")
# Ii from "Polymerase II" back to upper case
data$GeneSet <- str_replace_all(data$GeneSet, "Reactome", "Reactome:")
data$GeneSet <- str_replace_all(data$GeneSet, "Kegg", "KEGG:")
# GeneSet as factor
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)

# color palette
colfunc <- colorRampPalette(c("#D45E60", "#0F5057"))

# dotplot enrichment
g <- ggplot(data) +
  geom_point(aes(y = reorder(GeneSet,`[-log10FDR]`), 
                 x = `[-log10FDR]`, 
                 color = `Genes ratio`, 
                 size = `OR[log2(a*d)-log2(b*c)]`)) +
  geom_vline(xintercept = -log10(0.05),linetype="dashed", color = "black")+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_color_gradient(name = "% Genes ratio",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "OddsRatio") +
  xlab("-log10 FDR") +
  ylab("GOterm") +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
ggsave(path= "~/ownCloud/DexStim_RNAseq_Mouse/scripts/07_PlotsManuscript/plots_v3", "03_Figure3_enrichmentAbcd1_b.pdf", g, width = 8, height = 8)
