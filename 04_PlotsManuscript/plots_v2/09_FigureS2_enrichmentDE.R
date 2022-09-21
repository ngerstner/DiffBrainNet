##################################################
## Project: DexStim Mouse Brain
## Date: 02.03.2022
## Author: Nathalie
##################################################
# GO enrichment for DE genes across all regions

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)
library(writexl)


basepath <- "~/Documents/ownCloud/DexStim_RNAseq_Mouse/"


# 0. Read genes DE in all regions and background -----------------
genes <- read.table(file = paste0(basepath, "tables/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_entrezID.txt"),
                    header = FALSE)[,1]

# background are all genes in out dataset
background <- read.table(file = paste0(basepath, "tables/06_background_entrezID.txt"),
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


# 1.2 GO enrichment for upregulated genes DE in all regions ---------------------
genes_up <- read.table(file = paste0(basepath, "tables/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_up_entrezID.txt"),
                       header = FALSE)[,1]

# GO enrichment
min_genes <- round(length(genes_up)/10)
ego_up <- enrichGO(gene          = as.character(genes_up),
                   universe      = as.character(background),
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = min_genes,    # min number of genes associated with GO term
                   maxGSSize = 10000, # max number of genes associated with GO term
                   readable      = TRUE)@result
head(ego_up, n = 20)




# 1.3 GO enrichment for downregulated genes DE in all regions ---------------------
genes_down <- read.table(
  file = paste0(basepath, 
                "tables/06_overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_down_entrezID.txt"),
  header = FALSE)[,1]

# GO enrichment
min_genes <- length(genes_down)/10
ego_down <- enrichGO(gene          = as.character(genes_down),
                     universe      = as.character(background),
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = min_genes,    # min number of genes associated with GO term
                     maxGSSize = 10000, # max number of genes associated with GO term
                     readable      = TRUE)@result
head(ego_down, n = 20)


# 1.4 Merge dataframes from all, up and downregulated genes --------------------

data_go <- left_join(ego, ego_down, by = c("ID", "Description"), 
                     suffix = c(".all", ".down"))
data_go <- left_join(data_go, ego_up, by = c("ID", "Description"),
                     suffix = c("", ".up"))
# apparently not all entrez ids are valid --> not all taken into account in enrichment analysis
ngenes_rel <- 165
nback_rel <- 12089
# add columns relevant for gene ratio and odds ratio
data_go <- data_go %>%
  mutate(n_genes = as.numeric(str_extract(data_go$BgRatio.all, "^\\d+"))) %>%
  mutate(gr = Count.all/n_genes*100) %>%
  mutate(b = ngenes_rel - Count.all,
         c = n_genes - Count.all,
         d = nback_rel - Count.all,
         or = log2(Count.all*d)-log2(b*c))

data_heat <- data_go[1:25,c("Description", "p.adjust.down", "p.adjust")] %>%
  tidyr::pivot_longer(cols = p.adjust.down:p.adjust)
data_heat$Description <- factor(data_heat$Description,
                                levels = levels(reorder(data_go[1:25,]$Description, data_go[1:25,]$gr)))

# 1.5 Plot -------------------------------

#data_go[1:25,]$Description[order(data_go[1:25,]$Description)]

g <- ggplot(data_go[1:25,]) +
  geom_point(aes(y = reorder(Description,gr), 
                 #y = Description,
                 x = gr, 
                 color = -log10(p.adjust.all), 
                 size = or)) +
  scale_color_gradient(name = "-log10(FDR)",
                       low = "#D45E60", high = "#0F5057") +
  scale_size_continuous(name = "OddsRatio") +
  xlab("% GeneRatio") +
  ylab("GOterm") +
  theme_light() +
  theme(text = element_text(size= 14))

hm.1 <- 
  ggplot(data = data_heat, aes(
    x = Description,
    y = name,
    fill = value <= 0.01
  )) +
  geom_tile() +
  scale_y_discrete(name ="sig. GO term", 
                   limits=c("p.adjust","p.adjust.down"),
                   labels=c("upreg.", "downreg.")) +
  scale_fill_manual(
    name = "GO term sig.",
    values = c("darkgrey", "black") 
  ) +
  ylab("sig. GO term") +
  theme_light() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size =14),
    axis.text.y = element_text(size = 12),
    #legend.position = "top",
    legend.text = element_text(size = 10)
  ) +
  coord_flip() 

# combined plot
comb <- ggarrange(g, 
                  hm.1 +
                    theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank() ), 
                  nrow = 1,
                  widths = c(2.5,1))

# Save plot
ggexport(
  comb,
  filename = "09_FigureS2_enrichmentDE.pdf",
  height = 10,
  width = 13
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
tf$`[-log10 FDR]` <- -log10(tf$p.adjust.all)

# subset to columns that should be printed
tf <- tf %>% 
  dplyr::select(Category, GeneSet, n_genes, Count, b, c, d, gr, or, pvalue.all,
                p.adjust.all, `[-log10 FDR]`, geneID.all) %>%
  dplyr::rename(N_genes = n_genes,
                N_overlap_A = Count,
                `Input that does not overlap_B` = b,
                `GO term genes - overlap_C` = c,
                `Not GO term genes and not in my geneset_D` = d,
                `Genes ratio` = gr,
                `OR[log2(a*d)-log2(b*c)]` = or, 
                p = pvalue.all,
                adjP = p.adjust.all,
                genes = geneID.all)

# save to a excel file
write_xlsx(tf, "~/Documents/ownCloud/DexStim_RNAseq_Mouse/manuscript/Figures/Figure II_PFC/GOterms_DEgenesAllRegions.xlsx")
