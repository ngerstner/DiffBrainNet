##################################################
## Project: DexStim Mouse Brain
## Date: 30.11.2021
## Author: Nathalie
##################################################
# Expression level plots Syn3 and Abcd1 in PFC

library(stringr)
library(ggplot2)

reg <- "PFC"

ensembl_syn3 <- "ENSMUSG00000059602"
ensembl_abcd1 <- "ENSMUSG00000031378"

# read vsd normalized data 
data_exp <- read.table(paste0("~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_",
                              reg, "_deseq2_expression_vsd.txt"),
                       header = TRUE, sep = "\t")


### Syn3 ----------------

# subset Syn3 row
syn3 <- data_exp[ensembl_syn3,]

# get column indices for cntrl samples
indices_cntrl <- str_detect(colnames(data_exp), "CNTRL")

# df for cntrl samples
exp_cntrl <- data.frame("expression" = t(syn3[,indices_cntrl]),
                        "sample" = rownames(t(syn3[,indices_cntrl]))) %>%
  rename(expression = all_of(ensembl_syn3)) %>%
  mutate("treatment" = "CNTRL")

# df for dex samples
exp_dex <- data.frame("expression" = t(syn3[,!indices_cntrl]),
                      "sample" = rownames(t(syn3[,!indices_cntrl]))) %>%
  rename(expression = all_of(ensembl_syn3)) %>%
  mutate("treatment" = "DEX")

# combine cntrl and dex df
df <- rbind(exp_cntrl, exp_dex)
df$treatment <- factor(df$treatment)

ggplot(df, aes(x = treatment, y = expression, fill = treatment)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual("",
                    breaks = c("CNTRL", "DEX"),
                    labels = c("Baseline", "Treatment"),
                    values = c("#B0BFBB", "#46866E")) +
  scale_x_discrete(breaks = c("CNTRL", "DEX"),
                   labels = c("Baseline", "Treatment")) +
  theme_light() +
  theme(text = element_text(size= 16)) +
  xlab("") +
  ylab("Norm. Expression Level") +
  guides(fill = "none")

ggsave(filename = "09_FigureII_D.pdf", width = 7, height = 6)



### Abcd1 ----------------

# subset Abcd1 row
abcd1 <- data_exp[ensembl_abcd1,]

# get column indices for cntrl samples
indices_cntrl <- str_detect(colnames(data_exp), "CNTRL")

# df for cntrl samples
exp_cntrl <- data.frame("expression" = t(abcd1[,indices_cntrl]),
                        "sample" = rownames(t(abcd1[,indices_cntrl]))) %>%
  rename(expression = all_of(ensembl_abcd1)) %>%
  mutate("treatment" = "CNTRL")

# df for dex samples
exp_dex <- data.frame("expression" = t(abcd1[,!indices_cntrl]),
                      "sample" = rownames(t(abcd1[,!indices_cntrl]))) %>%
  rename(expression = all_of(ensembl_abcd1)) %>%
  mutate("treatment" = "DEX")

# combine cntrl and dex df
df <- rbind(exp_cntrl, exp_dex)
df$treatment <- factor(df$treatment)

ggplot(df, aes(x = treatment, y = expression, fill = treatment)) +
  geom_boxplot() +
  geom_jitter() +
  scale_fill_manual("",
                    breaks = c("CNTRL", "DEX"),
                    labels = c("Baseline", "Treatment"),
                    values = c("#B0BFBB", "#46866E")) +
  scale_x_discrete(breaks = c("CNTRL", "DEX"),
                   labels = c("Baseline", "Treatment")) +
  theme_light() +
  theme(text = element_text(size= 16)) +
  xlab("") +
  ylab("Norm. Expression Level") +
  guides(fill = "none")

ggsave(filename = "09_FigureII_E.pdf", width = 7, height = 6)
