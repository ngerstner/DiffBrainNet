##################################################
## Project: DexStim Mouse Brain
## Date: 29.11.2021
## Author: Nathalie
##################################################
# Pathway enrichment plot for manuscript 
# Tcf4 Differential Expression

library(stringr)
library(ggplot2)

regions <- c("AMY", "CER", "dCA1", "dDG", "PFC", "PVN", "vCA1", "vDG")

ensembl_tcf4 <- "ENSMUSG00000053477"

# Panel A: Expression values
df <- data.frame("region" = character(),
                 "treatment" = character(),
                 "sample" = character(),
                 "expression" = numeric())

for (reg in regions){
  
  # read vsd normalized data --> better raw data?
  data_exp <- read.table(paste0("~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
                                reg, "_deseq2_expression_vsd.txt"),
                        header = TRUE, sep = "\t")
  # subset Tcf4 row
  tcf4 <- data_exp[ensembl_tcf4,]
  
  # get column indices for cntrl samples
  indices_cntrl <- str_detect(colnames(data_exp), "CNTRL")
  
  # df for cntrl samples
  exp_cntrl <- data.frame("expression" = t(tcf4[,indices_cntrl]),
                          "sample" = rownames(t(tcf4[,indices_cntrl]))) %>%
    rename(expression = ENSMUSG00000053477) %>%
    mutate("treatment" = "CNTRL", "region" = reg)
  df <- rbind(df, exp_cntrl)
  
  # df for dex samples
  exp_dex <- data.frame("expression" = t(tcf4[,!indices_cntrl]),
                        "sample" = rownames(t(tcf4[,!indices_cntrl]))) %>%
    rename(expression = ENSMUSG00000053477) %>%
    mutate("treatment" = "DEX", "region" = reg)
  df <- rbind(df, exp_dex)
  
}

df$treatment <- factor(df$treatment)
df$region <- factor(df$region)

ggplot(df, aes(x = region, y = expression, fill = treatment)) +
  geom_boxplot() +
  scale_fill_manual("",
                    breaks = c("CNTRL", "DEX"),
                    labels = c("Baseline", "Treatment"),
                    values = c("#B0BFBB", "#46866E")) +
  theme_light() +
  theme(text = element_text(size= 14)) +
  xlab("Brain Region") +
  ylab("Norm. Expression Level")

ggsave(filename = "07_FigureSIV_A.pdf", width = 12, height = 8)


# Panel B: Fold Change in AMY, vDG and dDG

list_tcf4 <- list()

for (reg in c("AMY", "dDG", "vDG")){
# for (reg in regions){
  
  # read DE results
  data_exp <- read.table(paste0("~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
                                reg, "_deseq2_Dex_1_vs_0_lfcShrink.txt"),
                         header = TRUE, sep = "\t")
  tcf4 <- as.data.frame(data_exp[data_exp$Ensembl_ID == ensembl_tcf4,])
  print(tcf4)
  
  list_tcf4[[reg]] <- tcf4
  
}

df_tcf4 <- bind_rows(list_tcf4, .id = "region")


# bar plot Fold Change
fc <- ggplot(df_tcf4, aes(x = region, y = log2FoldChange) ) +
  geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity", fill = "darkgrey") +
  theme_light() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(text = element_text(size= 14)) +
  ylab("log2(FoldChange)") +
  xlab("Brain Region")

# barplot fdr p-value
fdr <- ggplot(df_tcf4, aes(x = region, y = -log10(padj)) ) +
  geom_bar(position =  position_dodge2(reverse=TRUE), stat="identity", fill = "#D45E60") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") + 
  theme_light() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(text = element_text(size= 14)) +
  ylab("-log10(FDR)") +
  xlab("")

# combined barplot
comb <- ggarrange(fc, 
                  fdr + 
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank() ), 
                  nrow = 1,
                  widths = c(1.3,1))

# ggexport(comb, filename = "07_FigureSIV_B_allRegions.pdf", width = 12, height = 8)
ggexport(comb, filename = "07_FigureSIV_B.pdf", width = 12, height = 8)

