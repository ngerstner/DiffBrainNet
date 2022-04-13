##################################################
## Project: DexStim Mouse Brain
## Date: 30.11.2021
## Author: Nathalie
##################################################
# Number of DE genes and unique percentage

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)

regions <-
  c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1")


# 1. Read DE tables from all regions ----------

list_reg_sig <- list()
list_genes_sig <- list()

for (reg in regions) {
  res <-
    fread(
      file = paste0(
        "~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables",
        "/02_",
        reg,
        "_deseq2_Dex_1_vs_0_lfcShrink.txt"
      ),
      sep = "\t"
    )
  na_indices <- which(is.na(res$padj))
  res$padj[na_indices] <- 1
  res_sig <- res[res$padj <= 0.1, ]
  # res_sig <- res[res$log2FoldChange >= 1]
  list_reg_sig[[reg]] <- res_sig
  list_genes_sig[[reg]] <- rownames(res_sig)
}

# 2. Concatenate DE tables -----------------

data <- bind_rows(list_reg_sig, .id = "region") %>%
  group_by(Ensembl_ID) %>%
  summarise(region = list(region))

data_unique <- data %>%
  mutate(nr_regions = lengths(region)) %>%
  mutate(unique = (nr_regions == 1)) %>%
  unnest(cols = c(region)) %>%
  mutate("combined_id" = paste0(region, "-", Ensembl_ID))

data_barplot <- data_unique %>%
  group_by(region, unique) %>%
  count() %>%
  group_by(region) %>%
  mutate(sum = sum(n))


# 3. Stacked barplot -------------------------

ggplot(data_barplot, aes(x = region, y = n, alpha = unique)) +
  geom_bar(position = "stack", stat = "identity", fill = "orange") +
  scale_alpha_manual(
    name = "",
    labels = c("DE in multiple regions", "DE unique"),
    values = c(1.0, 0.5)
  ) +
  xlab("Brain region") +
  ylab("# DE genes") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_blank(),
    legend.position = "top"
  ) +
  geom_text(aes(label = paste0(round((n / sum) * 100, digits = 1
  ), "%")),
  position = position_stack(vjust = 0.5),
  size = 4,
  show.legend = FALSE)
ggsave(
  "10_FigureSII_A.pdf",
  width = 8,
  height = 6
)
