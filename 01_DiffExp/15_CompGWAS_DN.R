##################################################
## Project: DexStim Mouse Brain
## Date: 07.09.2022
## Author: Nathalie
##################################################
# Compare DEs and DNs with previously identified GWAS risk genes 
# of 5 psychiatric disorders (https://pubmed.ncbi.nlm.nih.gov/32152537/)

# TODO: think about background --> restrict to genes in our dataset?

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(fgsea)
library(biomaRt)
library(ggplot2)

# trait of interest
traits <- c("ADHD", "ASD", "BD", "MDD", "SCZ")
# trait <- "MDD"

# define pathes
basedir <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse"
folder_table <- paste0(basedir,"/tables")
file_hmagma <- paste0(basedir, "/data/reviews/H-MAGMAv1.08_Adult_brain_output.xlsx")
files_dn <- list.files(path = paste0(folder_table, "/coExpression_kimono/03_AnalysisFuncoup/"), 
                       pattern = "_funcoup_treatment_nodebetweennessNorm_betacutoff0.01.csv$", full.names = TRUE)
regions_files <- sub(".*03a_(\\w*)_funcoup_treatment_nodebetween.*","\\1",files_dn)
file_background <- paste0(folder_table, "/06_background.txt")


## 1. Read DE genes from different brain regions
expr_list <- lapply(files_dn, function(x) fread(x) %>% filter(nodebetweenness_norm >= 1))
names(expr_list) <- regions_files

dn_list <- lapply(expr_list, function(x) x$ensembl_id)

# unique DN genes
data_unique <- bind_rows(expr_list, .id = "region") %>%
  group_by(ensembl_id) %>%
  summarise(region = list(region)) %>%
  mutate(nr_regions = lengths(region)) %>%
  mutate(unique = (nr_regions == 1)) %>%
  filter(unique)

dn_uni <- list()
for (r in regions_files){
  dn_uni[[r]] <- data_unique$ensembl_id[data_unique$region == r]
}

# background genes
background <- read.table(file_background, 
                         row.names = NULL)

# map mouse Ensembl Ids to human Ensembl Ids
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# --> this is a workaround as the normal commands lead to errors
# --> one mouse ID can be mapped to multiple human IDs

dn_human <- lapply(dn_list, function(x) getLDS(attributes=c("ensembl_gene_id"),
                                               filters="ensembl_gene_id", values=x, mart=mouse,
                                               attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1)

dn_uni_human <- lapply(dn_uni, function(x) getLDS(attributes=c("ensembl_gene_id"),
                                                  filters="ensembl_gene_id", values=x, mart=mouse,
                                                  attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1)

background_human <- getLDS(attributes=c("ensembl_gene_id"),
                           filters="ensembl_gene_id", values=background$V1, mart=mouse,
                           attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1


res_list <- list()
res_list_uni <- list()
for (trait in traits){

  ## 2. Read GWAS/H-MAGMA data
  data <- read_excel(file_hmagma, sheet=paste0("H-MAGMA_",trait)) %>%
    arrange(ZSTAT)  %>% # rank gene list according to zscore as suggested in paper
    filter(GENE %in% background_human)
  ranks <- data$ZSTAT
  names(ranks) <- data$GENE
  
  
  ## 3. Run fgsea for Gene Set Enrichment Analysis
  fgseaRes <- fgsea(pathways = dn_human, 
                    stats    = ranks) 
                    # sampleSize = 600)
  res_list[[trait]] <- fgseaRes
  
  # unique DN genes
  fgseaRes_uni <- fgsea(pathways = dn_uni_human, 
                        stats    = ranks)
  # sampleSize = 600)
  res_list_uni[[trait]] <- fgseaRes_uni
  
  # plotEnrichment(dn_human[["PFC"]],
  #                ranks) 

}

## 4. Plot GSEA p-values as heatmap from all traits and brain regions
df <- bind_rows(res_list, .id="trait")
df$FDR <- p.adjust(df$pval, method = "fdr") # apply correction across all tests
ggplot(df, aes(x = trait, y = pathway, fill = -log10(FDR))) +
  geom_tile() +
  xlab("GWAS trait") +
  ylab("Brain region") +
  scale_fill_gradient(name = "-log10(FDR)",
                      low="lightgrey", high="#D45E60",
                      limits = c(0,1)) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/15_Reviewer2_1b_DNall.pdf"),
       width = 8,
       height = 6)


# plot for unique DN genes
df <- bind_rows(res_list_uni, .id="trait")
df$FDR <- p.adjust(df$pval, method = "fdr") # apply correction across all tests
ggplot(df, aes(x = trait, y = pathway, fill = -log10(FDR))) +
  geom_tile() +
  xlab("GWAS trait") +
  ylab("Brain region") +
  scale_fill_gradient(name = "-log10(FDR)",
                      low="lightgrey", high="#D45E60",
                      limits = c(0,1)) +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12)
  )
ggsave(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/15_Reviewer2_1b_DNunique.pdf"),
       width = 8,
       height = 6)
