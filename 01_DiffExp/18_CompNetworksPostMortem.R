##################################################
## Project: DexStim Mouse Brain
## Date: 16.09.2022
## Author: Nathalie
##################################################
# Compare hub genes with disease associated modules
# from postmortem study
# (https://pubmed.ncbi.nlm.nih.gov/30545856/)


library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(biomaRt)
library(ggplot2)
library(ggpubr)


# define pathes
basedir <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse"
folder_table <- paste0(basedir,"/tables")
file_modules <- paste0(basedir, "/data/reviews/PostmortemBrain/aat8127_table_s5.xlsx")
file_de <- paste0(folder_table, "/02_PFC_deseq2_Dex_0_vs_1_lfcShrink.txt")
file_dn <- paste0(folder_table, "/coExpression_kimono/03_AnalysisFuncoup/", 
                  "03a_PFC_funcoup_treatment_nodebetweennessNorm_betacutoff0.01.csv")
file_background <- paste0(folder_table, "/06_background.txt")


## 1. Set up BioMart
# map mouse Ensembl Ids to human Ensembl Ids --> homology mapping
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# --> this is a workaround as the normal commands lead to errors

## 2. Read modules and disease associations
modules_genes <- read_excel(file_modules, sheet="geneModules")
modules_disease <- read_excel(file_modules, sheet = "Module-DiseaseAssociation") %>%
  filter(Network == "gene")


## 3. Read DE and hub genes from mouse in PFC
# DE genes
res <- read.table(file_de, sep="\t",
                  header = TRUE)
na_indices <- which(is.na(res$padj))
res$padj[na_indices] <- 1
res_sig <- res[res$padj <= 0.1,]
de_human <- getLDS(attributes=c("ensembl_gene_id"),
                    filters="ensembl_gene_id", values=rownames(res_sig), mart=mouse,
                    attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1
# hub genes
res_hub <- fread(file_dn) %>% filter(nodebetweenness_norm >= 1)
dn_human <- getLDS(attributes=c("ensembl_gene_id"),
                  filters="ensembl_gene_id", values=res_hub$ensembl_id, mart=mouse,
                  attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1

# DE and hub genes
dedn_human <- intersect(de_human, dn_human)
de_human <- setdiff(de_human, dedn_human)
dn_human <- setdiff(dn_human, dedn_human)

## 4. Read background mouse genes
background <- read.table(file_background, 
                         row.names = NULL)
background_human <- getLDS(attributes=c("ensembl_gene_id"),
                           filters="ensembl_gene_id", values=background$V1, mart=mouse,
                           attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1

## 5. DE and hub genes per module
modules_genes <- modules_genes %>%
  filter(ensembl_gene_id %in% background_human) %>% # subset modules to background genes
  mutate(de = ifelse(ensembl_gene_id %in% de_human, TRUE, FALSE), # add column to indicate if gene is DE gene
         hub = ifelse(ensembl_gene_id %in% dn_human, TRUE, FALSE),
         de_hub = ifelse(ensembl_gene_id %in% dedn_human, TRUE, FALSE)) # add column to indicate if gene is hub gene

# group by module and count
modules_count <- modules_genes %>%
  count(Module)

# count DE and hub per module
de_count <- modules_genes %>%
  group_by(Module) %>%
  count(de) %>%
  filter(de)

hub_count <- modules_genes %>%
  group_by(Module) %>%
  count(hub) %>%
  filter(hub)

dedn_count <- modules_genes %>%
  group_by(Module) %>%
  count(de_hub) %>%
  filter(de_hub)

# fisher's exact test for overrepresentation of DE and hub genes
# one-sided fisher's exact test with alternative "greater" is the same as 
# hypergeometric test (https://stats.stackexchange.com/questions/288081/use-fishers-exact-test-or-a-hypergeometric-test)
# Hypergeometric test (phyper) is also called by fora function from fgsea package
fis_test <- function(a1,b1,c1,d1){
  fisher.test(matrix(c(a1,b1,c1,d1), 
        nrow = 2, ncol = 2), alternative = "greater")$p.value
}
de_test <- modules_count %>%
  left_join(de_count[c("Module", "n")], by="Module") %>%
  replace(is.na(.), 0) %>%
  mutate(b = n.x - n.y,
         c = length(de_human) - n.y,
         d = nrow(modules_genes) - n.y - b - c) %>%
  rowwise() %>%
  mutate(fis_p = fis_test(n.y,b,c,d)) 
de_test$fis_fdr <- p.adjust(de_test$fis_p, method = "fdr")
de_test <- de_test %>%
  mutate(
    label_de = case_when(
      fis_fdr > 0.05 ~ "",
      fis_fdr > 0.01 ~ "*",
      fis_fdr > 0.001 ~ "**",
      !is.na(fis_fdr) ~ "***",
      TRUE ~ NA_character_
    )
  ) 

hub_test <- modules_count %>%
  left_join(hub_count[c("Module", "n")], by="Module") %>%
  replace(is.na(.), 0) %>%
  mutate(b = n.x - n.y,
         c = length(dn_human) - n.y,
         d = nrow(modules_genes) - n.y - b - c) %>%
  rowwise() %>%
  mutate(fis_p = fis_test(n.y,b,c,d)) 
hub_test$fis_fdr <- p.adjust(hub_test$fis_p, method = "fdr")
hub_test <- hub_test %>%
  mutate(
    label_hub = case_when(
      fis_fdr > 0.05 ~ "",
      fis_fdr > 0.01 ~ "*",
      fis_fdr > 0.001 ~ "**",
      !is.na(fis_fdr) ~ "***",
      TRUE ~ NA_character_
    )
  ) 

# join DE and hub counts to general count
df_count <- modules_count %>%
  left_join(de_count[c("Module", "n")], by="Module") %>%
  left_join(hub_count[c("Module", "n")], by="Module") %>%
  left_join(dedn_count[c("Module", "n")], by="Module") %>%
  mutate(o = factor(as.numeric(str_replace(Module, "geneM", "")),
                    levels = sort(as.numeric(str_replace(Module, "geneM", ""))))
         ) %>%
  replace(is.na(.), 0) %>%
  mutate(genes_de = n.y,
         genes_hub = n.x.x,
         genes_dehub = n.y.y,
         genes_rest = n.x-n.y-n.x.x-n.y.y) %>%
  pivot_longer(cols = starts_with("genes_"),
               names_to = "colour",
               values_to = "count") %>%
  mutate(colour = factor(colour, 
                         levels = c("genes_rest", "genes_hub", "genes_de", "genes_dehub"))) %>%
  left_join(de_test[c("Module", "label_de")], by="Module") %>%
  left_join(hub_test[c("Module", "label_hub")], by="Module")
df_count$label_de[df_count$colour != "genes_de"] <- ""
df_count$label_hub[df_count$colour != "genes_hub"] <- ""


# prepare disease df
df_disease <- modules_disease %>%
  mutate(o = factor(as.numeric(str_replace(ModulName, "geneM", "")),
                  levels = levels(df_count$o)))


## 6. Barplot
b <- ggplot(df_count, aes(x=o, y=count, fill=colour)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(breaks = c("genes_rest", "genes_hub", "genes_de", "genes_dehub"),
                    values = c("grey", "darkred", "orange", "#7E2274"),
                    labels = c("Module genes not DE or hub",
                               "Hub genes", "DE genes", "DE and hub genes")) +
  guides(fill = guide_legend(title = "Gene set")) +
  xlab("Module") +
  ylab("Number of genes") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

b_perc <- ggplot(df_count, aes(x=o, y=count, fill=colour)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label=label_de),
            position=position_fill(vjust=0.5), colour="white", size =6) +
  geom_text(aes(label=label_hub),
            position=position_fill(vjust=0.5), colour="white", size =6) +
  scale_fill_manual(breaks = c("genes_rest", "genes_hub", "genes_de", "genes_dehub"),
                    values = c("grey", "darkred", "orange", "#7E2274"),
                    labels = c("Module genes not DE or hub",
                               "Hub genes", "DE genes", "DE and hub genes")) +
  xlab("Module") +
  ylab("Percentage of genes") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

h <- ggplot(df_disease, aes(x=o, y=Group, fill=-log10(fdr))) +
  geom_tile() +
  scale_fill_gradient(name = "Disease association: -log10(FDR)",
                      low="lightgrey", high="#D45E60",
                      limits = c(0,max(-log10(df_disease$fdr)))) +
  scale_x_discrete(drop=FALSE) +
  ylab("Disease") +
  xlab("Module") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )


comb <- ggarrange(b + theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.x = element_blank()), 
                  b_perc + theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.title.x = element_blank()),
                  h , 
                  #legend.grob = get_legend(b),
                  nrow = 3,
                  heights=c(2,2,1),
                  align = "v")
comb

ggexport(comb, 
         filename = paste0(basedir, 
                           "/scripts_manuscript/04_PlotsManuscript/Revision/18_Reviewer2_2b_Networks_Percentage_significanceFDR.pdf"), 
         width = 14, height = 8)


