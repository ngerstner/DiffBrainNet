##################################################
## Project: DexStim Mouse Brain
## Date: 08.09.2020
## Author: Nathalie
##################################################
# Compare deseq SV DE genes between regions

setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

library(rlist)
library(VennDiagram)
library(UpSetR)
library(RColorBrewer)
library(org.Mm.eg.db)

# regions <- c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1", "vCA3", "dCA3")
# vCA3 and dCA3 were excluded
regions <- c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1")

folder_plots <- paste0("figures")
folder_tables <- paste0("tables")


# 0. functions -------------------------------
write_genelist <- function(genelist, filename){
  # write list with ENSEMBL IDs
  write.table(genelist, file = paste0(folder_tables, "/06_",filename,".txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  # write list with ENTREZ IDs
  entrez <- mapIds(org.Mm.eg.db, keys = genelist, keytype = "ENSEMBL", column="ENTREZID")
  write.table(entrez, file = paste0(folder_tables, "/06_",filename,"_entrezID.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  # write list with GENE SYMBOLS
  symbol <- mapIds(org.Mm.eg.db, keys = genelist, keytype = "ENSEMBL", column="SYMBOL")
  write.table(symbol, file = paste0(folder_tables, "/06_",filename,"_symbolID.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# 1. read DE tables from all regions ----------

list_reg <- list()
list_reg_sig <- list()
list_genes_sig <- list()
list_genes_sig_up <- list()
list_genes_sig_down <- list()
background_genes <- list()

for (reg in regions){
  res <- read.table(file=paste0(folder_tables, "/02_", reg, "_deseq2_Dex_1_vs_0_lfcShrink.txt"),sep="\t",
                    header = TRUE)
  na_indices <- which(is.na(res$padj))
  res$padj[na_indices] <- 1
  list_reg[[reg]] <- res
  res_sig <- res[res$padj <= 0.1,]
  list_reg_sig[[reg]] <- res_sig
  list_genes_sig[[reg]] <- res_sig$Ensembl_ID
  list_genes_sig_up[[reg]] <- res_sig$Ensembl_ID[res_sig$log2FoldChange > 0]
  list_genes_sig_down[[reg]] <- res_sig$Ensembl_ID[res_sig$log2FoldChange < 0]
  background_genes[[reg]] <- res$Ensembl_ID
}

background <- Reduce(intersect, background_genes)
write_genelist(background,"background")

union_back <- Reduce(union, background_genes)


# 2. Overlap between different results --------------

# png(filename = paste0(folder_plots, "/06_comparison_deseq_upsetPlot_withCA3.png"), height = 800, width = 1000)
png(filename = paste0(folder_plots, "/06_comparison_deseq_upsetPlot.png"), height = 600, width = 1000)
print(upset(fromList(list_genes_sig), nsets = 8, order.by = "freq",
      text.scale = c(1.8, 1.9, 1.8, 1.9, 1.9, 1.9),
      sets.x.label = "#DE genes in brain region",
      mainbar.y.label = "#DE genes in intersection"))
dev.off()

pdf(file = paste0(folder_plots, "/06_comparison_deseq_upsetPlot.pdf"), height = 9, width = 14)
print(upset(fromList(list_genes_sig), nsets = 8, order.by = "freq",
            text.scale = c(1.8, 1.9, 1.8, 1.9, 1.9, 1.9),
            sets.x.label = "#DE genes in brain region",
            mainbar.y.label = "#DE genes in intersection"))
dev.off()

# # 3. Overlap between regions --------------------
# # Regions AMY, CER, PFC and PVN
# venn.diagram(list(list_genes_sig[["AMY"]], list_genes_sig[["CER"]], list_genes_sig[["PFC"]], list_genes_sig[["PVN"]]),
#              category.names = c("AMY", "CER", "PFC", "PVN"),
#              filename = paste0(folder_plots, "/06_comparison_deseq_AMY_CER_PFC_PVN.png"),
#              output = TRUE,
#              imagetype="png" ,
#              height = 800, 
#              width = 800,
#              lwd = 1,
#              col=c("#284E5C","#288577","#73BA70","#EAE362"),
#              fill = c(alpha("#284E5C",0.3), alpha("#288577",0.3), alpha("#73BA70",0.3), alpha("#EAE362", 0.3)),
#              cex = 0.5,
#              fontfamily = "sans",
#              cat.cex = 0.3,
#              cat.default.pos = "outer",
#              cat.fontfamily = "sans",
#              cat.col = c("#284E5C","#288577","#73BA70","#EAE362"),
#              margin = 0.05)
# 
# # Hippocampus ventral
# venn.diagram(list(list_genes_sig[["vDG"]], list_genes_sig[["vCA1"]],
#                   list_genes_sig[["vCA3"]]),
#              category.names = c("vDG", "vCA1", "vCA3"),
#              filename = paste0(folder_plots, "/06_comparison_deseq_HIPventral.png"),
#              output = TRUE,
#              imagetype="png" ,
#              height = 800, 
#              width = 800,
#             lwd = 1,
#             col=c("#284E5C","#288577","#73BA70"),
#             fill = c(alpha("#284E5C",0.3), alpha("#288577",0.3), alpha("#73BA70",0.3)),
#             cex = 0.5,
#             fontfamily = "sans",
#             cat.cex = 0.3,
#             cat.default.pos = "outer",
#             cat.fontfamily = "sans",
#             cat.col = c("#284E5C","#288577","#73BA70"),
#             margin = 0.05)
# 
# # Hippocampus dorsal
# venn.diagram(list(list_genes_sig[["dDG"]], list_genes_sig[["dCA1"]],
#                   list_genes_sig[["dCA3"]]),
#              category.names = c("dDG", "dCA1", "dCA3"),
#              filename = paste0(folder_plots, "/06_comparison_deseq_HIPdorsal.png"),
#              output = TRUE,
#              imagetype="png" ,
#              height = 800, 
#              width = 800,
#              lwd = 1,
#              col=c("#284E5C","#288577","#73BA70"),
#              fill = c(alpha("#284E5C",0.3), alpha("#288577",0.3), alpha("#73BA70",0.3)),
#              cex = 0.5,
#              fontfamily = "sans",
#              cat.cex = 0.3,
#              cat.default.pos = "outer",
#              cat.fontfamily = "sans",
#              cat.col = c("#284E5C","#288577","#73BA70"),
#              margin = 0.05)

# PFC, PVN and dCA1 for progress report on 26.11.2020
library(eulerr)
png(filename = paste0(folder_plots, "/06_comparison_deseq_dCA1_PFC_PVN.png"),
    height = 600, width = 600)
plot(euler(list("dCA1" = list_genes_sig[["dCA1"]], "PFC" = list_genes_sig[["PFC"]],
                "PVN" = list_genes_sig[["PVN"]]), shape = "ellipse"), 
     labels = list(cex = 1.5), quantities = list(cex = 1.5))
dev.off()

  
# 4. Gene lists -------------------------------------

# 4.1 Overlap all regions 
overlap <- Reduce(intersect, list_genes_sig)
write_genelist(overlap, "overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1")
overlap_up <- Reduce(intersect, list_genes_sig_up)
write_genelist(overlap_up, "overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_up")
overlap_down <- Reduce(intersect, list_genes_sig_down)
write_genelist(overlap_down, "overlap_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1_down")

union_sig <- Reduce(union, list_genes_sig)
write_genelist(union_sig, "union_AMY-CER-PFC-PVN-dDG-vDG-dCA1-vCA1")

# 4.2 Overlap all regions except PFC
overlap1 <- Reduce(intersect, list_genes_sig[c(1,3:8)])
write_genelist(overlap1, "overlap_AMY-CER-PVN-dDG-vDG-dCA1-vCA1")

# difference between all regions and all regions without PFC
d <- setdiff(overlap1, overlap)
d <- list_reg[["PFC"]][d,]
d$ENTREZID <- mapIds(org.Mm.eg.db, keys = rownames(d), keytype = "ENSEMBL", column = "ENTREZID")
d$SYMBOL <- mapIds(org.Mm.eg.db, keys = rownames(d), keytype = "ENSEMBL", column = "SYMBOL")
write.table(d, file = paste0(folder_tables, "/06_diff_allRegionsExceptPFC.txt"),
            quote = FALSE)

# 4.3 Overlap of PFC and CER
overlap_pfc_cer <- Reduce(intersect, list_genes_sig[c(2,4)])
union_other <- Reduce(union, list_genes_sig[c(1,3,5:8)])
write_genelist(overlap_pfc_cer, "overlap_CER-PFC")
d <- setdiff(overlap_pfc_cer, union_other)
write_genelist(d, "overlap_CER-PFC_only")


# 4.3 Overlap of dorsal hippocampus
overlap_dhip <- Reduce(intersect, list_genes_sig[c(6,8)])
overlap_dhip_up <- Reduce(intersect, list_genes_sig_up[c(6,8)])
overlap_dhip_down <- Reduce(intersect, list_genes_sig_down[c(6,8)])
write_genelist(overlap_dhip, "overlap_dDG-dCA1")


# 4.4 Overlap of ventral hippocampus
overlap_vhip <- Reduce(intersect, list_genes_sig[c(5,7)])
overlap_vhip_up <- Reduce(intersect, list_genes_sig_up[c(5,7)])
overlap_vhip_down <- Reduce(intersect, list_genes_sig_down[c(5,7)])
write_genelist(overlap_vhip, "overlap_vDG-vCA1")


# 4.5 overlap ventral and dorsal HIP  -------------
venn.diagram(
  list(overlap_dhip, overlap_vhip),
  category.names = c("dorsal HIP", "ventral HIP"),
  filename = paste0(folder_plots, "/06_comparison_deseq_HIP.png"),
  output = TRUE,
  imagetype = "png" ,
  height = 800,
  width = 800,
  lwd = 1,
  col = c("#284E5C", "#288577"),
  fill = c("#284E5C", "#288577"),
  alpha = 0.3,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.dist = 0.05,
  cat.fontfamily = "sans",
  cat.col = c("#284E5C", "#288577"),
  margin = 0.05
)

venn.diagram(
  list(overlap_dhip_up, overlap_vhip_up, overlap_dhip_down, overlap_vhip_down),
  category.names = c("dorsal HIP up", "ventral HIP up", "dorsal HIP down", "ventral HIP down"),
  filename = paste0(folder_plots, "/06_comparison_deseq_HIP_upDown.png"),
  output = TRUE,
  imagetype = "png" ,
  height = 800,
  width = 800,
  lwd = 1,
  col=c("#284E5C","#288577","#73BA70","#EAE362"),
  fill = c("#284E5C","#288577","#73BA70","#EAE362"),
  alpha = 0.3,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  # cat.default.pos = "text",
  cat.dist = 0.1,
  cat.fontfamily = "sans",
  cat.col = c("#284E5C", "#288577", "#73BA70", "#EAE362"),
  margin = 0.05)

# intersection dorsal ventral
int <- intersect(overlap_dhip, overlap_vhip)
write_genelist(int, "HIP_intersectDorsalVentral")
# DE in dorsal but not ventral
diff_dv <- setdiff(overlap_dhip, overlap_vhip)
write_genelist(diff_dv, "HIP_diffDorsalVentral")
# DE in ventral but not dorsal
diff_vd <- setdiff(overlap_vhip, overlap_dhip)
write_genelist(diff_vd, "HIP_diffVentralDorsal")
# union dorsal ventral
un_dv <- union(overlap_dhip, overlap_vhip)
write_genelist(un_dv, "HIP_unionDorsalVentral")
