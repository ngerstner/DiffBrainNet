##################################################
## Project: DexStim Mouse Brain
## Date: 29.11.2021
## Author: Nathalie
##################################################
# GWAS enrichment plot for manuscript 
# vCA1

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

# 1. Pathway enrichment for Abcd1 and neighbours

data <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure III_vCA1/network vCA1/Unique DEs and their diff neighbours_FUMA_gene2func67313/vCA1 Unique DE genes+diff neighbours_FUMA enrichments.xlsx",
                  sheet = "GWAS")

data <- arrange(data, desc(`[-log10 of FDR`))
data <- data[data$`my adj p BH` <= 0.05,]
data <- data[!is.na(data$`my adj p BH`),]

# GeneSet as factor 
data$GeneSet <- factor(data$GeneSet, levels = data$GeneSet)

# bar plot gene ratio
gr <- ggplot(data, aes(x = GeneSet, y = `Genes ratio`) ) +
  geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity", fill = "#0F5057") +
  scale_fill_manual() +
  theme_light() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(text = element_text(size= 14)) +
  ylab("% Gene Ratio") +
  labs(x = NULL)

# bar plot odds ratio
or <- ggplot(data, aes(x = GeneSet, y = `OR[log2(a*d)-log2(b*c)]`) ) +
  geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity", fill = "#FAA916") +
  scale_fill_manual() +
  theme_light() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(text = element_text(size= 14)) +
  ylab("Odds Ratio") +
  labs(x = NULL)

# barplot fdr p-value
fdr <- ggplot(data, aes(x = GeneSet, y = `[-log10 of FDR`) ) +
  geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity", fill = "#D45E60") +
  geom_hline(yintercept = -log10(0.05),linetype="dashed", color = "red") + 
  theme_light() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(text = element_text(size= 14)) +
  ylab("-log10(FDR)") +
  labs(x = NULL)

# combined barplot
comb <- ggarrange(gr, 
                  or + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank() ),
                  fdr + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank() ), 
                  nrow = 1,
                  widths = c(1.8,1,1))

ggexport(comb, filename = "06_FigureSIII_A.pdf", width = 14, height = 8)
