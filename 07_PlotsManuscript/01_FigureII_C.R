##################################################
## Project: DexStim Mouse Brain
## Date: 21.11.2021
## Author: Nathalie
##################################################
# GO Enrichment plot for manuscript 
# Unique DE and hub genes in PFC

# set working directory to source file location
setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse/scripts/07_PlotsManuscript")

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)


# 1. GO enrichment unique DE and hub genes in PFC

data_hub <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/FUMA_diff unique hubgenes_gene2func67664/PFC diff unique hubgenes_FUMA.xlsx",
                     sheet = "GO bp") 

data_DE <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/DE PFC/FUMA_gene2func66314/PFC unique DE genes _FUMA.xlsx",
                      sheet = "GO bp") 


# data for left panel with top 14 terms for DE
top_DE <- data_DE %>%
  dplyr::filter(N_overlap_A > 24)
top_DE <- top_DE[order(top_DE$`my adj p BH`, decreasing = FALSE),][1:14,]

top_DE_hub <- data_hub %>%
  filter(GeneSet %in% top_DE$GeneSet)

top_DE_comb <- bind_rows("de" = top_DE, "hub" = top_DE_hub, .id = "group")


# data for right panel with top 14 terms for hubs
top_hub <- data_hub %>%
  dplyr::filter(N_overlap_A > 2)
top_hub <- top_hub[order(top_hub$`my adj p BH`, decreasing = FALSE),][1:14,]

top_hub_DE <- data_DE %>%
  filter(GeneSet %in% top_hub$GeneSet)

missing_term <- setdiff(top_hub$GeneSet, top_hub_DE$GeneSet)
add_entry <- data.frame("GeneSet" = missing_term)
top_hub_DE <- bind_rows(top_hub_DE, add_entry)

top_hub_comb <- bind_rows("hub" = top_hub, "de" = top_hub_DE, .id = "group")


# combine everything into one df for plotting
top_data <- bind_rows("de" = top_DE_comb, "hub" = top_hub_comb, .id = "panel")
# remove the "GO_" from GO terms
top_data$GeneSet <- str_replace_all(top_data$GeneSet, '^GO_', ' ')
# remove the "_" from GO terms
top_data$GeneSet <- str_replace_all(top_data$GeneSet, "_", " ")
# only capitalize first letter per word
top_data$GeneSet <- str_to_title(top_data$GeneSet)
# of and to should start with lower case
top_data$GeneSet <- str_replace_all(top_data$GeneSet, "Of", "of")
top_data$GeneSet <- str_replace_all(top_data$GeneSet, "To", "to")

# labeller function for facet titles
facet_names <- c(
  'de'="Top 14 GO terms for DE genes",
  'hub'="Top 14 GO terms for hub genes"
)

# barplot
g <- ggplot(top_data, aes(x = GeneSet, y = `[-log10 fdr`, fill = group)) +
        geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity") +
         
        scale_colour_manual(values = c("grey", "black"), guide = FALSE) +
        scale_fill_manual(name = "Geneset",
                          values = c("orange", "darkred"),
                          labels = c("DE genes", "hub genes")) +
        coord_flip() +
        scale_x_discrete(limits = rev) +
        xlab("GOterm") +
        ylab("-log10(FDR)") +
        ggtitle("GO terms enriched for unique DE and hub genes in PFC") +
        facet_wrap(~panel, scales="free", labeller = as_labeller(facet_names)) +
        theme_light() +
        theme(text = element_text(size= 14))
print(g)
ggsave("01_FigureII_C.pdf", g, width = 14, height = 8)
  
