prior <- read.table(
  file = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/prior_expr_funcoup_mm.csv",
  sep = ",",
  header = TRUE)

background <- read.table(
  file = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/06_background.txt"
)

de_PFC <- read.table(
  file = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/02_AMY_deseq2_Dex_1_vs_0_lfcShrink.txt",
  header = TRUE
) %>%
  filter(padj <= 0.1)


prior_genes <- unique(c(prior$V1, prior$V2))

intersect(prior_genes, background$V1)

intersect(prior_genes, de_PFC$Ensembl_ID)

prior <- prior %>%
  filter(Gene_A %in% background$V1 & Gene_B %in% background$V1)


funcoup <- read.table(
  file = "~/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/FC5.0_M.musculus_compact",
  sep = "\t"
)
