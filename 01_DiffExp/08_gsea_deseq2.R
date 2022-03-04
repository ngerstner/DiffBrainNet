##################################################
## Project: DexStim Mouse Brain
## Date: 08.09.2020
## Author: Nathalie
##################################################
# Gene Set Enrichment Analysis of genes with high fold change
# !!! too many gene sets tested --> nothing remains sig after multiple testing !!!

library(fgsea)
library(org.Mm.eg.db)

region <- "PFC"

folder_plots <- paste0("figures")
folder_tables <- paste0("tables")

res <- read.table(file=paste0(folder_tables, "/02_", region, "_deseq2_Dex_0_vs_1_lfcShrink_lfc0.5.txt"),sep="\t")
res <- res[order(res$log2FoldChange),]
ranks <- res$log2FoldChange
entrez <- mapIds(org.Mm.eg.db, keys = rownames(res), keytype = "ENSEMBL", column="ENTREZID")
names(ranks) <- entrez
head(ranks)


barplot(sort(ranks, decreasing = T))
data(examplePathways)

fgseaRes <- fgsea(examplePathways, ranks)

head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes <- fgseaRes[order(padj, -abs(NES)),]

topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-pval)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-pval)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(examplePathways[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

