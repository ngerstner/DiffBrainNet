##################################################
## Project: DexStim Mouse Brain
## Date: 13.04.2021
## Author: Nathalie
##################################################
# Functional annotation with clusterProfiler
# make figure for manuscript

library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(enrichplot)
library(gridExtra)
library(stringr)


basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"



# 0. Read genes DE in all regions and background -----------------
genes <- read.table(file = paste0(basepath, "tables/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_entrezID.txt"),
                    header = FALSE)[,1]

# background are all genes in out dataset
background <- read.table(file = paste0(basepath, "tables/06_background_entrezID.txt"),
                         header = FALSE)[,1]



# 1.1 GO enrichment for genes DE in all regions ---------------------

# GO enrichment
# TODO: decide on maxGSSize --> with 10000 very similar results to anRichment
# --> Anthi and me decided that it makes sense to leave the cutoff very high
# (no point of restricting the terms here)
ego <- enrichGO(gene          = as.character(genes),
                universe      = as.character(background),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                minGSSize = 10,    # min number of genes associated with GO term
                maxGSSize = 10000, # max number of genes associated with GO term
                readable      = TRUE)
head(ego, n = 20)

barplot(ego, showCategory=20)
dotplot(ego, showCategory=30) + ggtitle("dotplot for DE genes in all regions")

# SIMPLIFY enriched GO terms (remove very similar terms)
ego_simple <- clusterProfiler::simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)@result
head(ego_simple, n = 20)

#barplot(ego_simple, showCategory=20)
#dotplot(ego_simple, showCategory=30) + ggtitle("dotplot for DE genes in all regions")



# 1.2 GO enrichment for upregulated genes DE in all regions ---------------------
genes_up <- read.table(file = paste0(basepath, "tables/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_up_entrezID.txt"),
                    header = FALSE)[,1]

# GO enrichment
# TODO: decide on maxGSSize --> with 10000 very similar results to anRichment
ego_up <- enrichGO(gene          = as.character(genes_up),
                universe      = as.character(background),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                minGSSize = 10,    # min number of genes associated with GO term
                maxGSSize = 10000, # max number of genes associated with GO term
                readable      = TRUE)@result
head(ego_up, n = 20)




# 1.3 GO enrichment for downregulated genes DE in all regions ---------------------
genes_down <- read.table(
  file = paste0(basepath, 
                "tables/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_down_entrezID.txt"),
  header = FALSE)[,1]

# GO enrichment
# TODO: decide on maxGSSize --> with 10000 very similar results to anRichment
ego_down <- enrichGO(gene          = as.character(genes_down),
                   universe      = as.character(background),
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   minGSSize = 10,    # min number of genes associated with GO term
                   maxGSSize = 10000, # max number of genes associated with GO term
                   readable      = TRUE)@result
head(ego_down, n = 20)


# 1.4 Merge dataframes from all, up and downregulated genes --------------------

data_go <- left_join(ego_simple, ego_down, by = c("ID", "Description"), 
                     suffix = c(".all", ".down"))
data_go <- left_join(data_go, ego_up, by = c("ID", "Description"),
                     suffix = c("", ".up"))

data_heat <- data_go[1:30,c("Description", "p.adjust.down", "p.adjust")] %>%
  tidyr::pivot_longer(cols = p.adjust.down:p.adjust)


# 1.5 Barplot -------------------------------

bp.1 <-
  ggplot(data = data_go[1:25,], aes(
    x = factor(Description, levels = rev(data_go$Description[1:25])),
    y = Count.all/165
  )) +
  geom_bar(stat = "identity", position = "stack",
           fill = "#226666") +
  # scale_fill_manual(
  #   name = "",
  #   labels = c("DE in multiple regions", "DE unique"),
  #   values = c("red3", "navy")
  # ) +
  scale_y_continuous(trans="reverse") +
  ylab("Gene ration") +
  xlab("GO terms - biological process") +
  theme_light() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size =14),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 10)
  ) +
  coord_flip() 

bp.2 <- 
  ggplot(data = data_go[1:25,], aes(
    x = factor(Description, levels = rev(data_go$Description[1:25])),
    y = -log10(p.adjust.all)
  )) +
  geom_bar(stat = "identity", position = "stack",
           fill = "#AA3939") +
  # scale_fill_manual(
  #   name = "",
  #   labels = c("DE in multiple regions", "DE unique"),
  #   values = c("red3", "navy")
  # ) +
  ylab("-log10(adj. p-value)") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 10)
  ) +
  coord_flip() 

hm.1 <- 
  ggplot(data = data_heat, aes(
    x = factor(Description, levels = rev(data_go$Description[1:25])),
    y = name,
    fill = value <= 0.01
    # fill = p.adjust
  )) +
  geom_tile() +
  scale_y_discrete(name ="sig. GO term", 
                   limits=c("p.adjust","p.adjust.down"),
                  labels=c("upreg.", "downreg.")) +
  scale_fill_manual(
    name = "GO term significant",
    values = c("darkgrey", "#FFB620") 
  ) +
  ylab("sig. GO term") +
  theme_light() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    # legend.position = "none"
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  coord_flip() 

plot_comb <- grid.arrange(bp.1, bp.2, hm.1, nrow = 1,
                          widths = c(3, 1, 1.5))

# Save plot
ggsave(
  plot_comb,
  filename = paste0(basepath, "figures/07a_goEnrichment_allRegions.png"),
  height = 10,
  width = 12
)

# 2.1 Disease gene enrichment for genes DE in all regions ---------------------

# Map ENTREZ IDs from mouse to human
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 <- getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", 
                  values = genes , mart = mouse, attributesL = c("entrezgene_id"), 
                  martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])

backgroundV2 <- getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", 
                     values = background , mart = mouse, attributesL = c("entrezgene_id"), 
                     martL = human, uniqueRows=T)
humanb <- unique(backgroundV2[,2])

# Disease enrichment
# TODO: decide on maxGSSize
dgn <- enrichDGN(gene          = as.character(humanx),
                universe      = as.character(humanb),
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                minGSSize = 500,    # min number of genes associated with GO term
                maxGSSize = 5000, # max number of genes associated with GO term
                readable      = TRUE)
head(dgn@result$Description, n = 100)

dgn_x <- pairwise_termsim(dgn)
enrichplot::emapplot(dgn_x, showCategory = 50)

dgn@result[dgn@result$Description == "Schizophrenia",]


# x <- enrichDO(gene          = as.character(humanx),
#               ont           = "DO",
#               pvalueCutoff  = 0.05,
#               pAdjustMethod = "BH",
#               universe      = as.character(humanb),
#               minGSSize     = 5,
#               maxGSSize     = 500,
#               qvalueCutoff  = 0.05,
#               readable      = FALSE)
# head(x)
# 
# x2 <- pairwise_termsim(dgn)
# enrichplot::emapplot(x2, showCategory = 50)


library(disgenet2r)
library(psygenet2r)

data2 <- gene2disease(
  gene     = humanx,
  vocabulary = "ENTREZ",
  # database = "PSYGENET",
  score =c(0.2, 1),
  verbose  = TRUE
)

data2_table <- data2@qresult
data2_table <- data2_table[(data2_table$disease_class_name == "   Mental Disorders" |
                             data2_table$disease_class_name == "   Nervous System Diseases"),]
data2_table <- data2_table[(str_detect(data2_table$disease_class_name, "Mental Disorders") |
                              str_detect(data2_table$disease_class_name, "Nervous System Diseases")),]

plot( data2,
      class = "Network",
      prop = 10)

plot( data2,
      class  ="Heatmap")

plot( data2,
      class="DiseaseClass")

# disease enrichment using disgenet2r
# does not work for PSYGENET as database (bug in code)
enr <- disease_enrichment(
  genes = humanx,
  universe = humanb,
  vocabulary = "ENTREZ",
  verbose = TRUE,
  database = "CURATED",
  warnings = TRUE
)@qresult

# gene-disease associations (GDA) using psygenet2r
m1 <- psygenetGene(
  gene     = humanx, 
  database = "ALL",
  verbose  = TRUE
)
plot( m1, type = "GDCA network" )
plot( m1 )
plot( m1, type="GDCA heatmap" )
# geneAttrPlot( m1, type = "disease category", class = "Lollipop" )

png(filename = paste0(basepath, "figures/07a_diseaseAssociations_allRegions.png"),
    width = 600, height = 600)
plot( m1 )
dev.off()

# disease enrichment (per disease class) using psygenet2r
enr_psy <- enrichedPD(
  gene = humanx,
  verbose = TRUE,
  warnings = TRUE
)

ggplot(enr_psy, aes(x = MPD,
                    y = -log10(p.value))) +
  geom_bar(stat = "identity",
           position = "stack",
           fill = "#1E88E5") +
  coord_flip()
ggsave(
  filename = paste0(basepath, "figures/07a_diseaseEnrichment_allRegions.png"),
  width = 8,
  height = 8
)
