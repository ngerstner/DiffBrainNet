##################################################
## Project: DexStim Mouse Brain
## Date: 08.09.2020
## Author: Nathalie
##################################################
# Functional annotation with anRichment

library(anRichment)
library(anRichmentMethods)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)

basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"
folder_plots <- paste0("figures")
folder_tables <- paste0("tables")

GOcollection <- buildGOcollection(organism = "mouse")

# 1.1 GO enrichment all regions ---------------------
genes <- read.table(file = paste0(basepath, folder_tables, "/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_entrezID.txt"),
                  header = FALSE)
background <- read.table(file = paste0(basepath, folder_tables, "/06_background_entrezID.txt"),
                         header = FALSE)
modules <- rep("not_significant", nrow(background))
modules[which(background$V1 %in% genes$V1)] <- "significant"

# enrichment
GOenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = background$V1,
  refCollection = GOcollection,
  useBackground = "given",
  # threshold = 0.1,
  # thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "not_significant"
)

# Filter out terms with less than 10 genes overlap and FDR above 0.1
table.display <- GOenrichment$enrichmentTable %>%
  filter(nCommonGenes >= 10 & FDR <= 0.1)
write.csv(table.display, file = paste0(folder_tables,"/07_GOenrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.csv"),
          row.names = FALSE)

# Plot results (Top 20 terms)
table.display %>%
  head(20) %>%
  arrange(desc(FDR)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
ggplot(aes(x=dataSetName, y = -log10(FDR), fill = inGroups)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  scale_fill_discrete(name = "Ontology",
                      labels = c("BP", "MF", "CC")) +
  coord_flip() +
  ylab("-log10(FDR)") +
  xlab("GOterm") +
  ggtitle("GO terms enriched for DE genes across all brain regions (Top 20)") +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 15, hjust=1))
ggsave(filename = paste0(folder_plots, "/07_GOenrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.png"))


# 1.2 NCBI BioSystems collection -----------------
biosysCollection <- BioSystemsCollection("mouse")
knownGroups(biosysCollection)

# enrichment
KEGGenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = background$V1,
  refCollection = biosysCollection,
  useBackground = "given",
  threshold = 0.1,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "not_significant"
)
collectGarbage()

names(KEGGenrichment)
names(KEGGenrichment$enrichmentTable)
table.display <- KEGGenrichment$enrichmentTable
table.display$overlapGenes <- shortenStrings(table.display$overlapGenes, maxLength = 70,
                                             split = "|")
write.csv(GOenrichment$enrichmentTable, file = paste0(folder_tables,"/07_PathwayEnrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.csv"),
          row.names = FALSE)

# plot results (plotted pvalues are not adjusted for multiple testing)
table.display %>%
  arrange(desc(pValue)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
ggplot(aes(x=dataSetName, y = -log10(pValue), fill=inGroups)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  coord_flip() +
  ggtitle("Pathways enriched for DE genes across all brain regions")
ggsave(filename = paste0(folder_plots, "/07_PathwayEnrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.png"))


# 1.3 calculate enrichment of disease associated genes ---------
library(disgenet2r)
# Schizophrenia C0036341
# Major depressive disorder C1269683
# Bipolar disorder C0005586
genes_SCZ <- disease2gene( disease = c("C0036341"), 
                         database = "CURATED", verbose = TRUE )@qresult
genes_MDD <- disease2gene( disease = c("C1269683"), 
                           database = "CURATED", verbose = TRUE )@qresult
genes_BIP <- disease2gene( disease = c("C0005586"), 
                           database = "CURATED", verbose = TRUE )@qresult

genesymbols <- read.table(file = paste0(folder_tables, "/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_symbolID.txt"),
                    header = FALSE)
genesymbols <- sapply(genesymbols$V1, toupper)

bgsymbols <- read.table(file = paste0(folder_tables, "/06_background_symbolID.txt"),
                          header = FALSE)
bgsymbols <- sapply(bgsymbols$V1, toupper)

pval <- data.frame(disease = c("SCZ", "MDD", "BIP"), pvalue = c(1,1,1))
diseaseGenes <- data.frame(disease = character(), genes = character())
# hypergeometric test SCZ
tmp <- intersect(genes_SCZ$gene_symbol, genesymbols)
diseaseGenes <- rbind(diseaseGenes, cbind(rep("SCZ", times = length(tmp)), tmp))
o <- length(intersect(genes_SCZ$gene_symbol, genesymbols))
o_b <- length(intersect(genes_SCZ$gene_symbol, bgsymbols))
test <- fisher.test(matrix(c(o, o_b-o, length(genesymbols)-o, length(bgsymbols)-length(genesymbols)-o_b+o), 2, 2),
                    alternative='greater')$p.value
pval[1,2] <- test

# hypergeometric test MDD
tmp <- intersect(genes_MDD$gene_symbol, genesymbols)
diseaseGenes <- rbind(diseaseGenes, cbind(rep("MDD", times = length(tmp)), tmp))
o <- length(intersect(genes_MDD$gene_symbol, genesymbols))
o_b <- length(intersect(genes_MDD$gene_symbol, bgsymbols))
test <- fisher.test(matrix(c(o, o_b-o, length(genesymbols)-o, length(bgsymbols)-length(genesymbols)-o_b+o), 2, 2),
                    alternative='greater')$p.value
pval[2,2] <- test

# hypergeometric test BIP
tmp <- intersect(genes_BIP$gene_symbol, genesymbols)
diseaseGenes <- rbind(diseaseGenes, cbind(rep("BIP", times = length(tmp)), tmp))
o <- length(intersect(genes_BIP$gene_symbol, genesymbols))
o_b <- length(intersect(genes_BIP$gene_symbol, bgsymbols))
test <- fisher.test(matrix(c(o, o_b-o, length(genesymbols)-o, length(bgsymbols)-length(genesymbols)-o_b+o), 2, 2),
                    alternative='greater')$p.value
pval[3,2] <- test

pval[,2] <- p.adjust(pval[,2], method = "fdr")

# plot disease gene enrichment
ggplot(pval, aes(x=disease, y=-log10(pvalue))) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  coord_flip() +
  xlab("disease") +
  ylab("-log10(FDR)") +
  ggtitle("Enrichment of known disease genes among DE genes shared between all regions") +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 15))
ggsave(filename = paste0(folder_plots, "/07_diseaseGeneEnrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.png"))

write.table(diseaseGenes, file = paste0(folder_tables, "/07_diseaseGenes_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(pval, file = paste0(folder_tables, "/07_diseaseGenes_pvalues_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1.txt"),
            quote = FALSE, row.names = FALSE)


# 2. Hippocampus: dorsal vs ventral ------------

intersect_hip <- read.table(file = paste0(folder_tables, "/06_HIP_intersectDorsalVentral_entrezID.txt"),
                            header = FALSE)
diff_dorsal <- read.table(file = paste0(folder_tables, "/06_HIP_diffDorsalVentral_entrezID.txt"),
                           header = FALSE)
diff_ventral <- read.table(file = paste0(folder_tables, "/06_HIP_diffVentralDorsal_entrezID.txt"),
                           header = FALSE)
union_hip <- read.table(file = paste0(folder_tables, "/06_HIP_unionDorsalVentral_entrezID.txt"),
                        header = FALSE)
# !!!! BACKGROUND are only diff exp genes in any HIP area here !!!!
# union_hip <- background
modules <- rep("XXX", nrow(union_hip))
modules[which(union_hip$V1 %in% intersect_hip$V1)] <- "intersect"
modules[which(union_hip$V1 %in% diff_dorsal$V1)] <- "dorsal"
modules[which(union_hip$V1 %in% diff_ventral$V1)] <- "ventral"

# enrichment
GOenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = union_hip$V1,
  refCollection = GOcollection,
  useBackground = "given",
  nBestDataSets = length(GOcollection$dataSets),
  # threshold = 0.1,
  # thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE
)

# Filter to the BP terms and exclude terms with less 10 genes overlap
# Exclude terms with nominal pvalue above 0.1
table.display <- GOenrichment$enrichmentTable %>%
  filter(inGroups == "GO|GO.BP|GO" & class != "XXX") %>% 
  filter(nCommonGenes >= 10 & pValue <= 0.1) %>%
  group_by(class) %>% slice_min(order_by = FDR, n = 10)
write.csv(GOenrichment$enrichmentTable, file = paste0(folder_tables,"/07_GOenrichment_HIP.csv"),
          row.names = FALSE)

# Plot results (plotted pvalues are not adjusted for multiple testing)
table.display %>%
  # arrange(desc(class)) %>%
  mutate(dataSetName=factor(dataSetName,levels=rev(dataSetName))) %>%
ggplot(aes(x=dataSetName, y = -log10(FDR), fill = class)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  coord_flip() +
  xlab("GOterm") +
  ggtitle("GO terms enriched for DE genes in Hippocampus")
ggsave(filename = paste0(folder_plots, "/07_GOenrichment_HIP.png"))


# 2.2 NCBI BioSystems collection -----------------
biosysCollection <- BioSystemsCollection("mouse")
knownGroups(biosysCollection)

# enrichment
KEGGenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = union_hip$V1,
  refCollection = biosysCollection,
  useBackground = "given",
  threshold = 0.1,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE
)
collectGarbage()

names(KEGGenrichment)
names(KEGGenrichment$enrichmentTable)
table.display <- KEGGenrichment$enrichmentTable
table.display$overlapGenes <- shortenStrings(table.display$overlapGenes, maxLength = 70,
                                             split = "|")
write.csv(GOenrichment$enrichmentTable, file = paste0(folder_tables,"/07_PathwayEnrichment_HIP.csv"),
          row.names = FALSE)

# plot results (plotted pvalues are not adjusted for multiple testing)
table.display %>%
  arrange(desc(class)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
ggplot(aes(x=dataSetName, y = -log10(pValue), fill=class)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  coord_flip() +
  ggtitle("Pathways enriched for DE genes in Hippocampus")
ggsave(filename = paste0(folder_plots, "/07_PathwayEnrichment_HIP.png"))




# 3.1 GO enrichment all regions without PFC ---------------------
# !!! not really interesting --> 9 genes from which 5 are nominally sig also in PFC !!!!
GOcollection <- buildGOcollection(organism = "mouse")
genes <- read.table(file = paste0(folder_tables, "/06_diff_allRegionsExceptPFC.txt"),
                    header = TRUE)
background <- read.table(file = paste0(folder_tables, "/06_background_entrezID.txt"),
                         header = FALSE)
modules <- rep("not_significant", nrow(background))
modules[which(background$V1 %in% genes$ENTREZID)] <- "significant"

# enrichment
GOenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = background$V1,
  refCollection = GOcollection,
  useBackground = "given",
  # threshold = 0.1,
  # thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "not_significant"
)

table.display <- GOenrichment$enrichmentTable
write.csv(table.display, file = paste0(folder_tables,"/07_GOenrichment_AMY-CER-PVN-dDG-vDG-dCA1-vCA1.csv"),
          row.names = FALSE)

# plot results (plotted pvalues are not adjusted for multiple testing)
table.display %>%
  arrange(desc(pValue)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
  ggplot(aes(x=dataSetName, y = -log10(pValue), fill = inGroups)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  scale_fill_discrete(name = "Ontology",
                      labels = c("BP", "MF", "CC")) +
  coord_flip() +
  ylab("-log10(pValue)") +
  xlab("GOterm") +
  ggtitle("GO terms enriched for DE genes across all brain regions except PFC") +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 15, hjust=1))
ggsave(filename = paste0(folder_plots, "/07_GOenrichment_AMY-CER-PVN-dDG-vDG-dCA1-vCA1.png"))




# 4.1 GO enrichment only PFC and CER ---------------------
genes <- read.table(file = paste0(folder_tables, "/06_overlap_CER-PFC_only_entrezID.txt"),
                    header = FALSE)
background <- read.table(file = paste0(folder_tables, "/06_background_entrezID.txt"),
                         header = FALSE)
modules <- rep("not_significant", nrow(background))
modules[which(background$V1 %in% genes$V1)] <- "significant"

# enrichment
GOenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = background$V1,
  refCollection = GOcollection,
  useBackground = "given",
  # threshold = 0.1,
  # thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "not_significant"
)

table.display <- GOenrichment$enrichmentTable
write.csv(table.display, file = paste0(folder_tables,"/07_GOenrichment_CER-PFC.csv"),
          row.names = FALSE)

# plot results (plotted pvalues are not adjusted for multiple testing)
table.display %>%
  arrange(desc(pValue)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
  ggplot(aes(x=dataSetName, y = -log10(pValue), fill = inGroups)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  scale_fill_discrete(name = "Ontology",
                      labels = c("BP", "MF", "CC")) +
  coord_flip() +
  ylab("-log10(pValue)") +
  xlab("GOterm") +
  ggtitle("GO terms enriched for DE genes across all brain regions except PFC") +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 15, hjust=1))
ggsave(filename = paste0(folder_plots, "/07_GOenrichment_AMY-CER-PVN-dDG-vDG-dCA1-vCA1.png"))




# 5.1 GO enrichment all regions upregulated---------------------
genes <- read.table(file = paste0(folder_tables, "/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_up_entrezID.txt"),
                    header = FALSE)
background <- read.table(file = paste0(folder_tables, "/06_background_entrezID.txt"),
                         header = FALSE)
modules <- rep("not_significant", nrow(background))
modules[which(background$V1 %in% genes$V1)] <- "significant"

# GO enrichment
GOenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = background$V1,
  refCollection = GOcollection,
  useBackground = "given",
  nBestDataSets = length(GOcollection$dataSets),
  # threshold = 0.1,
  # thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "not_significant"
)

# Exclude terms with less than 10 genes overlap or FDR above 0.1
table.display <- GOenrichment$enrichmentTable %>%
  filter(nCommonGenes >= 10 & FDR <= 0.1)
write.csv(table.display, file = paste0(folder_tables,"/07_GOenrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_up.csv"),
          row.names = FALSE)

# Plot results (Top 20)
table.display %>%
  head(20) %>%
  arrange(desc(FDR)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
  ggplot(aes(x=dataSetName, y = -log10(FDR), fill = inGroups)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  scale_fill_discrete(name = "Ontology",
                      labels = c("BP", "MF", "CC")) +
  coord_flip() +
  ylab("-log10(FDR)") +
  xlab("GOterm") +
  ggtitle("GO terms enriched for upregulated DE genes across all brain regions (Top 20)") +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 15, hjust=1))
ggsave(filename = paste0(folder_plots, "/07_GOenrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_up.png"))



# 5.2 GO enrichment all regions downregulated---------------------
genes <- read.table(file = paste0(folder_tables, "/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_down_entrezID.txt"),
                    header = FALSE)
background <- read.table(file = paste0(folder_tables, "/06_background_entrezID.txt"),
                         header = FALSE)
modules <- rep("not_significant", nrow(background))
modules[which(background$V1 %in% genes$V1)] <- "significant"

# GO enrichment
GOenrichment <- enrichmentAnalysis(
  classLabels = modules,
  identifiers = background$V1,
  refCollection = GOcollection,
  useBackground = "given",
  nBestDataSets = length(GOcollection$dataSets),
  # threshold = 0.1,
  # thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "not_significant"
)

# Exclude terms with less than 10 genes overlap or FDR above 0.1
table.display <- GOenrichment$enrichmentTable %>%
  filter(nCommonGenes >= 10 & FDR <= 0.1)
write.csv(table.display, file = paste0(folder_tables,"/07_GOenrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_down.csv"),
          row.names = FALSE)

# Plot results (Top 20)
table.display %>%
  head(20) %>%
  arrange(desc(FDR)) %>%
  mutate(dataSetName=factor(dataSetName,levels=dataSetName)) %>%
  ggplot(aes(x=dataSetName, y = -log10(FDR), fill = inGroups)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = -log10(0.1),linetype="dashed", color = "red") +
  scale_fill_discrete(name = "Ontology",
                      labels = c("BP", "MF", "CC")) +
  coord_flip() +
  ylab("-log10(FDR)") +
  xlab("GOterm") +
  ggtitle("GO terms enriched for downregulated DE genes across all brain regions (Top 20)") +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 15, hjust=1))
ggsave(filename = paste0(folder_plots, "/07_GOenrichment_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_down.png"))