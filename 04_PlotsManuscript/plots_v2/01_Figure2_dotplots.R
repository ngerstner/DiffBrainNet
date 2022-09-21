##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GO Enrichment plot for manuscript 
# Unique DE and hub genes in PFC

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)


# 1. GO enrichment unique DE and hub genes in PFC

data_hub <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/networks/FUMA_diff unique hubgenes_gene2func67664/PFC diff unique hubgenes_FUMA.xlsx",
                      sheet = "GO bp") 

data_DE <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/DE PFC/FUMA_gene2func66314/PFC unique DE genes _FUMA.xlsx",
                     sheet = "GO bp") 


# data for left panel with top 14 terms for DE
top_DE <- data_DE %>%
  dplyr::filter(N_overlap_A > 24)
top_DE <- top_DE[order(top_DE$`my adj p BH`, decreasing = FALSE),][1:14,]

top_DE_hub <- data_hub %>%
  filter(GeneSet %in% top_DE$GeneSet)

top_DE_comb <- bind_rows("DE genes" = top_DE, "hub genes" = top_DE_hub, .id = "group")


# data for right panel with top 14 terms for hubs
top_hub <- data_hub %>%
  dplyr::filter(N_overlap_A > 2)
top_hub <- top_hub[order(top_hub$`my adj p BH`, decreasing = FALSE),][1:14,]

top_hub_DE <- data_DE %>%
  filter(GeneSet %in% top_hub$GeneSet)

missing_term <- setdiff(top_hub$GeneSet, top_hub_DE$GeneSet)
add_entry <- data.frame("GeneSet" = missing_term)
top_hub_DE <- bind_rows(top_hub_DE, add_entry)

top_hub_comb <- bind_rows("hub genes" = top_hub, "DE genes" = top_hub_DE, .id = "group")


clean <- function(top_data) {
  # remove the "GO_" from GO terms
  top_data$GeneSet <- str_replace_all(top_data$GeneSet, '^GO_', ' ')
  # remove the "_" from GO terms
  top_data$GeneSet <- str_replace_all(top_data$GeneSet, "_", " ")
  # only capitalize first letter per word
  top_data$GeneSet <- str_to_title(top_data$GeneSet)
  # of and to should start with lower case
  top_data$GeneSet <- str_replace_all(top_data$GeneSet, "Of", "of")
  top_data$GeneSet <- str_replace_all(top_data$GeneSet, "To", "to")
  
  return(top_data)
}

top_DE_comb <- clean(top_DE_comb)
top_hub_comb <- clean(top_hub_comb)

# order of terms
top_DE_comb$GeneSet <- factor(top_DE_comb$GeneSet, levels = top_DE_comb$GeneSet[14:1])

# dotplot terms DE
g <- ggplot(top_DE_comb) +
  geom_point(aes(#y = reorder(GeneSet,-`my adj p BH`), 
                 y = GeneSet,
                 x = group, 
                 color = `[-log10 fdr`, 
                 size = `Genes ratio`)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "GeneRatio") +
  xlab("") +
  ylab("GOterm") +
  #ggtitle("GO enrichment for unique DE genes in PFC") +
  # facet_wrap(~panel, scales="free", labeller = as_labeller(facet_names)) +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
ggsave("01_Figure2_GO_DE.pdf", g, width = 8, height = 8)


# order of terms
top_hub_comb$GeneSet <- factor(top_hub_comb$GeneSet, levels = top_hub_comb$GeneSet[14:1])

# dotplot terms hub
g <- ggplot(top_hub_comb) +
  geom_point(aes(y = reorder(GeneSet,-`my adj p BH`), 
                 x = group, 
                 color = `[-log10 fdr`, 
                 size = `Genes ratio`)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "GeneRatio") +
  xlab("") +
  ylab("GOterm") +
  #ggtitle("GO enrichment for unique hub genes in PFC") +
  # facet_wrap(~panel, scales="free", labeller = as_labeller(facet_names)) +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
ggsave("01_Figure2_GO_hub.pdf", g, width = 8, height = 8)

