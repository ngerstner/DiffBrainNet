##################################################
## Project: DexStim Mouse Brain
## Date: 29.11.2021
## Author: Nathalie
##################################################
# GO enrichment plot for manuscript 
# Tcf4 baseline network in PFC

library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)

# 1. GO enrichment for vCA unique DE genes and their diff neighbours

data <- read_xlsx(path = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure IV_GWAS/TCF4/FUMA_baselin Tcf4PFCnetwork_gene2func68958/FUMA Tcf4BaselinePFCnetwork.xlsx",
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
fdr <- ggplot(data, aes(x = GeneSet, y = `[-log10]`) ) +
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
                  widths = c(1.9,1,1))

ggexport(comb, filename = "05_FigureIV_Bbase.pdf", width = 14, height = 10)


