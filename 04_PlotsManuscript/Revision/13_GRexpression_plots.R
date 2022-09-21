##################################################
## Project: DexStim Mouse Brain
## Date: 06.09.2022
## Author: Nathalie
##################################################
# Plot GR expression level at control level in all brain regions
# and check if correlation with number of DE genes exists
# Suggestion of reviewer

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

basedir <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse"
folder_table <- paste0(basedir,"/tables")
mean_exp_file <- file.path(folder_table, "12_meanExp_regionDex.csv")
files <- list.files(path = folder_table, pattern = "deseq2_Dex_0_vs_1_lfcShrink.txt$", full.names = TRUE)

# IDs
gene_symbol_gr <- "Nr3c1"
gene_symbol_mr <- "Nr3c2"
ensembl_gr <- "ENSMUSG00000024431"
ensembl_mr <- "ENSMUSG00000031618"


## 1. Read mean expression data
data <- read.table(mean_exp_file, sep = ";", header = TRUE) %>%
  dplyr::select(contains("CNTRL"), gene_symbol) %>%
  dplyr::rename_all(funs(str_replace_all(., "_CNTRL", "")))
gr_mr <- data[data$gene_symbol == gene_symbol_gr | data$gene_symbol == gene_symbol_mr,] %>% # subset to GR  and MR
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "gene_symbol")
gr_mr <- as.data.frame(t(gr_mr)) %>%
  tibble::rownames_to_column("region")

## 2. Number of DE genes per region
# read files with DE results
regions_files <- sub(".*02_(\\w*)_deseq2.*","\\1",files)
expr_list <- lapply(files, function(x) fread(x))
names(expr_list) <- regions_files
# concatenate DE results to get those for GR and MR
de_df <- bind_rows(expr_list, .id="region") %>% 
  dplyr::filter(V1 == ensembl_gr | V1 == ensembl_mr) %>%
  dplyr::mutate(gene = case_when(
    V1 == ensembl_gr ~ "Nr3c1",
    V1 == ensembl_mr ~ "Nr3c2"
  )) %>%
  dplyr::select(region, gene, log2FoldChange, padj) %>%
  pivot_wider(names_from = gene, values_from = c(log2FoldChange,padj))
# filter each df to genes with FDR smaller 0.1
expr_list <- lapply(expr_list, function(x) dplyr::filter(x, padj <= 0.1))
df_degenes <- bind_rows(expr_list, .id="region") %>%
  mutate(up_down = ifelse(log2FoldChange > 0, "up", "down")) 
n_degenes <- df_degenes %>%
  group_by(region) %>%
  count()
n_updowngenes <- df_degenes %>%
  group_by(region, up_down) %>%
  count()

## 3. Join mean expression values with number of DE genes
gr_mr <- inner_join(gr_mr, n_degenes, by = "region")
gr_mr <- dplyr::rename(gr_mr, nr_de = n)
gr_mr <- inner_join(gr_mr, de_df, by="region")

gr_mr_updown <- inner_join(gr_mr, n_updowngenes, by = "region") %>%
  mutate(up_down = factor(gr_mr_updown$up_down, levels = c("up", "down")))

## 4. Plot expression of GR and MR against number of DE genes
# GR
ggplot(gr_mr, aes(x=Nr3c1, y=nr_de, 
                  size=log2FoldChange_Nr3c1, 
                  color=-log10(padj_Nr3c1),
                  label=region)) +
  geom_point() +
  geom_text(hjust=-0.25, vjust=-0.25, 
            size=4, color="black") +
  xlim(min(gr_mr$Nr3c1),max(gr_mr$Nr3c1)+0.1) +
  xlab("Normalised Expression Level GR") +
  ylab("Number of DE genes") +
  scale_color_gradient(name = "-log10(FDR)",
                       low="lightgrey", high="#D45E60",
                       limits = c(0,13)) +
  scale_size_continuous(name="log2(FoldChange)") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/13_Reviewer1_GR.pdf"),
       width = 8,
       height = 6)

# MR
ggplot(gr_mr, aes(x=Nr3c2, y=nr_de, 
                  size=log2FoldChange_Nr3c2, 
                  color=-log10(padj_Nr3c2),
                  label=region)) +
  geom_point() +
  geom_text(hjust=-0.25, vjust=-0.25, 
            size=4, color="black") +
  xlim(min(gr_mr$Nr3c2),max(gr_mr$Nr3c2)+0.1) +
  xlab("Normalised Expression Level MR") +
  ylab("Number of DE genes") +
  scale_color_gradient(name = "-log10(FDR)",
                       low="lightgrey", high="#D45E60",
                       limits = c(0,13)) +
  scale_size_continuous(name="log2(FoldChange)") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/13_Reviewer1_MR.pdf"),
       width = 8,
       height = 6)


## 5. Plot expression of GR and MR against number of DE genes - stratified by up and downregulation
# GR
ggplot(gr_mr_updown, aes(x=Nr3c1, y=n, 
                  size=log2FoldChange_Nr3c1, 
                  color=-log10(padj_Nr3c1),
                  label=region)) +
  geom_point() +
  geom_text(hjust=-0.25, vjust=-0.25, 
            size=4, color="black") +
  xlim(min(gr_mr_updown$Nr3c1),max(gr_mr_updown$Nr3c1)+0.1) +
  xlab("Normalised Expression Level GR") +
  ylab("Number of DE genes") +
  scale_color_gradient(name = "-log10(FDR)",
                       low="lightgrey", high="#D45E60",
                       limits = c(0,13)) +
  scale_size_continuous(name="log2(FoldChange)") +
  facet_wrap(~up_down, nrow = 2, scales = "free",
             labeller = labeller(up_down = 
                                   c("down" = "Downregulated genes",
                                     "up" = "Upregulated genes")
             )) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/13_Reviewer1_GR_updown.pdf"),
       width = 8,
       height = 6)

# MR
ggplot(gr_mr_updown, aes(x=Nr3c2, y=n, 
                  size=log2FoldChange_Nr3c2, 
                  color=-log10(padj_Nr3c2),
                  label=region)) +
  geom_point() +
  geom_text(hjust=-0.25, vjust=-0.25, 
            size=4, color="black") +
  xlim(min(gr_mr_updown$Nr3c2),max(gr_mr_updown$Nr3c2)+0.1) +
  xlab("Normalised Expression Level MR") +
  ylab("Number of DE genes") +
  scale_color_gradient(name = "-log10(FDR)",
                       low="lightgrey", high="#D45E60",
                       limits = c(0,13)) +
  scale_size_continuous(name="log2(FoldChange)") +
  facet_wrap(~up_down, nrow = 2, scales = "free",
             labeller = labeller(up_down = 
                                   c("down" = "Downregulated genes",
                                     "up" = "Upregulated genes")
             )) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/13_Reviewer1_MR_updown.pdf"),
       width = 8,
       height = 6)


## 6. Plot GR/MR ratio against number of DE genes - stratified by up and downregulation
gr_mr_updown <- gr_mr_updown %>%
  mutate(grmr_ratio = Nr3c1/Nr3c2)

# GR/MR ratio
ggplot(gr_mr_updown, aes(x=grmr_ratio, y=n, 
                         size=log2FoldChange_Nr3c2, 
                         color=-log10(padj_Nr3c2),
                         label=region)) +
  geom_point() +
  geom_text(hjust=-0.25, vjust=-0.25, 
            size=4, color="black") +
  xlim(min(gr_mr_updown$grmr_ratio),max(gr_mr_updown$grmr_ratio)+0.1) +
  xlab("Normalised Expression Level MR") +
  ylab("Number of DE genes") +
  scale_color_gradient(name = "-log10(FDR)",
                       low="lightgrey", high="#D45E60",
                       limits = c(0,13)) +
  scale_size_continuous(name="log2(FoldChange)") +
  facet_wrap(~up_down, nrow = 2, scales = "free",
             labeller = labeller(up_down = 
                                   c("down" = "Downregulated genes",
                                     "up" = "Upregulated genes")
             )) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/13_Reviewer1_GRMRratio_updown.pdf"),
       width = 8,
       height = 6)
