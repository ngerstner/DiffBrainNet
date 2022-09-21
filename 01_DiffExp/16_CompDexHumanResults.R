##################################################
## Project: DexStim Mouse Brain
## Date: 14.09.2022
## Author: Nathalie
##################################################
# Compare DEs and DNs with DE genes from human Dex study
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8669026/)


library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(readxl)
library(biomaRt)
library(ggplot2)
library(UpSetR)
library(rlist)


# define pathes
basedir <- "/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse"
folder_table <- paste0(basedir,"/tables")
file_human <- paste0(basedir, "/data/reviews/DexHumanBlood_MooreEtAl/41398_2021_1756_MOESM1_ESM.xlsx")
files_de <- list.files(path = folder_table, pattern = "deseq2_Dex_0_vs_1_lfcShrink.txt$", full.names = TRUE)
print(sub(".*02_(\\w*)_deseq2.*","\\1",files_de))
regions_files <- sub(".*02_(\\w*)_deseq2.*","\\1",files_de)
file_background <- paste0(folder_table, "/06_background.txt")

## 0. Set up BioMart
# map mouse Ensembl Ids to human Ensembl Ids --> homology mapping
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# --> this is a workaround as the normal commands lead to errors


## 1. Read differentially expressed genes from human data into list
data <- read_excel(file_human, sheet="2_DEA_results") 

# get Ensembl IDs for human genes
ens_human <- getBM(attributes = c("ensembl_gene_id", "illumina_humanht_12_v3"), # no difference between v3 and v4
                   filters = "illumina_humanht_12_v3",
                   values = data$Probe_Id, mart = human)
data <- dplyr::inner_join(data, ens_human, by = c("Probe_Id"="illumina_humanht_12_v3"))

## 2.Read background genes
background <- read.table(file_background, 
                         row.names = NULL)
# get human ensembl IDs for background genes
background_human <- getLDS(attributes=c("ensembl_gene_id"),
                    filters="ensembl_gene_id", values=background$V1, mart=mouse,
                    attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1

## 3. Subset human DE genes
# subset human DE genes to those in background set
data <- data[data$ensembl_gene_id %in% background_human,]

# DE genes
data_de <- data[data$padj <= 0.05,]

## 4. Read DE tables from all regions 
expr_list <- lapply(files_de, function(x) fread(x) %>% filter(padj <= 0.1))
names(expr_list) <- regions_files

de_list <- lapply(expr_list, function(x) x$V1)
de_list_human <- lapply(de_list, 
                        function(x) getLDS(attributes=c("ensembl_gene_id"),
                                           filters="ensembl_gene_id", values=x, mart=mouse,
                                           attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1)

# unique DE genes
data_unique <- bind_rows(expr_list, .id = "region") %>%
  group_by(V1) %>%
  summarise(region = list(region)) %>%
  mutate(nr_regions = lengths(region)) %>%
  mutate(unique = (nr_regions == 1)) %>%
  filter(unique)

de_uni <- list()
for (r in regions_files){
  tmp <- data_unique$V1[data_unique$region == r]
  de_uni[[r]] <- getLDS(attributes=c("ensembl_gene_id"),
                        filters="ensembl_gene_id", values=tmp, mart=mouse,
                        attributesL=c("ensembl_gene_id"), martL=human)$Gene.stable.ID.1
}


## 5. Plot overlap between human and mouse DE genes
# unique DE genes
list_flat <- c(de_uni, list("human" = data_de$ensembl_gene_id))
pdf(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/16_Reviewer2_0_DEGDexHuman_DEunique.pdf"),
    width = 12, height = 8)
par(mar=c(6,5,4,2) + 0.1)
print(upset(fromList(list_flat), nsets = 9, order.by = "freq",
            text.scale = c(1.8, 1.9, 1.8, 1.9, 1.9, 1.9),
            sets.x.label = "#DE genes",
            mainbar.y.label = "#DE genes in intersection"))
dev.off()

# all DE genes
list_flat <- c(de_list_human, list("human" = data_de$ensembl_gene_id))
pdf(paste0(basedir, "/scripts_manuscript/04_PlotsManuscript/Revision/16_Reviewer2_0_DEGDexHuman_DEall.pdf"),
    width = 12, height = 8)
print(upset(fromList(list_flat), nsets = 9, order.by = "freq",
            text.scale = c(1.8, 1.9, 1.8, 1.9, 1.9, 1.9),
            sets.x.label = "#DE genes",
            mainbar.y.label = "#DE genes in intersection"))
dev.off()
