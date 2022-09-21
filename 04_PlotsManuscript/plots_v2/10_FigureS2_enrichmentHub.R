##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GO enrichment for hub genes across all regions


library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggpubr)
library(writexl)

regions <-
  c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1")


# 1. Read DE tables from all regions ----------

list_reg_sig <- list()
list_genes_sig <- list()

for (reg in regions) {
  res <-
    fread(
      file = paste0(
        "~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/coExpression_kimono/03_AnalysisFuncoup/",
        "/04_",
        reg,
        "_funcoup_differential_nodebetweennessNorm_betacutoff0.01.csv"
      ),
      sep = ","
    )
  na_indices <- which(is.na(res$nodebetweenness_norm))
  res$padj[na_indices] <- 0
  res_sig <- res[res$nodebetweenness_norm >= 1.0, ]
  list_reg_sig[[reg]] <- res_sig
  list_genes_sig[[reg]] <- res_sig$ensembl_id
}

hub_shared <- Reduce(intersect, list_genes_sig)
genes <- mapIds(org.Mm.eg.db, keys = hub_shared, keytype = "ENSEMBL", column="ENTREZID")

# background are all genes in out dataset
background <- read.table(file = paste0("~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/06_background_entrezID.txt"),
                         header = FALSE)[,1]


# 1.1 GO enrichment for genes DE in all regions ---------------------

# GO enrichment
min_genes <- round(length(genes)/10)
ego <- enrichGO(gene          = as.character(genes),
                universe      = as.character(background),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                minGSSize = min_genes,    # min number of genes associated with GO term
                maxGSSize = 10000, # max number of genes associated with GO term
                readable      = TRUE)@result
head(ego, n = 20)


# 1.4 Merge dataframes from all, up and downregulated genes --------------------

# apparently not all entrez ids are valid --> not all taken into account in enrichment analysis
ngenes_rel <- 7
nback_rel <- 12089
# add columns relevant for gene ratio and odds ratio
data_go <- ego %>%
  mutate(n_genes = as.numeric(str_extract(ego$BgRatio, "^\\d+"))) %>%
  mutate(gr = Count/n_genes*100) %>%
  mutate(b = ngenes_rel - Count,
         c = n_genes - Count,
         d = nback_rel - Count,
         or = log2(Count*d)-log2(b*c))

data_go$or[data_go$or == Inf] <-12

# 1.5 Plot -------------------------------

g <- ggplot(data_go[1:25,]) +
  geom_point(aes(y = reorder(Description,gr), 
    #y = Description,
    x = gr, 
    color = -log10(p.adjust), 
    size = or)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "OddsRatio") +
  xlab("% GeneRatio") +
  ylab("GOterm") +
  theme_light() +
  theme(text = element_text(size= 14))
print(g)
# Save plot
ggsave(
  "10_FigureS2_enrichmentHub.pdf",
  g,
  height = 10,
  width = 10
)


# bring table to a format similar as other tables
tf <- data_go[1:25,]

# add Category column
tf$Category <- "GO_bp"
# add GeneSet column
tf$GeneSet <- tf$Description %>% 
  stringr::str_replace_all(" ", "_") %>%
  stringr::str_to_upper() 
tf$GeneSet <- paste0("GO_", tf$GeneSet)
# add -log10 transformed FDR pvalue
tf$`[-log10 FDR]` <- -log10(tf$p.adjust)

# subset to columns that should be printed
tf <- tf %>% 
  dplyr::select(Category, GeneSet, n_genes, Count, b, c, d, gr, or, pvalue,
                p.adjust, `[-log10 FDR]`, geneID) %>%
  dplyr::rename(N_genes = n_genes,
                N_overlap_A = Count,
                `Input that does not overlap_B` = b,
                `GO term genes - overlap_C` = c,
                `Not GO term genes and not in my geneset_D` = d,
                `Genes ratio` = gr,
                `OR[log2(a*d)-log2(b*c)]` = or, 
                p = pvalue,
                adjP = p.adjust,
                genes = geneID)

# save to a excel file
write_xlsx(tf, "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/GOterms_hubsAllRegions.xlsx")

