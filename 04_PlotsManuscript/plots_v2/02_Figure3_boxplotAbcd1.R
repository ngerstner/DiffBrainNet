##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# Expression level plot Abcd1 in all brain regions


library(stringr)
library(ggplot2)

regions <- c("AMY", "CER", "dCA1", "dDG", "PFC", "PVN", "vCA1", "vDG")

ensembl_abcd1 <- "ENSMUSG00000031378"

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
  # subset Abcd1 row
  abcd1 <- data_exp[ensembl_abcd1,]
  
  # get column indices for cntrl samples
  indices_cntrl <- str_detect(colnames(data_exp), "CNTRL")
  
  # df for cntrl samples
  exp_cntrl <- data.frame("expression" = t(abcd1[,indices_cntrl]),
                          "sample" = rownames(t(abcd1[,indices_cntrl]))) %>%
    rename(expression = ENSMUSG00000031378) %>%
    mutate("treatment" = "CNTRL", "region" = reg)
  df <- rbind(df, exp_cntrl)
  
  # df for dex samples
  exp_dex <- data.frame("expression" = t(abcd1[,!indices_cntrl]),
                        "sample" = rownames(t(abcd1[,!indices_cntrl]))) %>%
    rename(expression = ENSMUSG00000031378) %>%
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

ggsave(filename = "02_Figure3_boxplotAbcd1.pdf", width = 12, height = 8)
