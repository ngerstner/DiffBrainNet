##################################################
## Project: DexStim Mouse Brain
## Date: 03.08.2020
## Author: Nathalie
##################################################
# RNA-seq Analysis. Part 2. Analysis with DESeq2
# Brain region AMYGDALA
# !!!!This is the right script, use this one!!!!

setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

# 1. load data and packages ----
source("scripts/01_DiffExp/01_setup_NG.R")
source("scripts/01_DiffExp/00_functions.R")
library(MASS)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lattice)
library(EnhancedVolcano)
library(magrittr)
library(PCAtools)
library(org.Mm.eg.db)

region <- "PFC"
prefix_plots <- paste0("figures/02_", region,"_deseq2_")
prefix_tables <- paste0("tables/02_",region,"_deseq2_")

# 1. Subset data to region of interest ----
Rcovariates <- droplevels(subset(covariates, Region==region))
Rcountdata <- countdata[,rownames(Rcovariates)]


# 2. Create DESeq object ----
all(Rcovariates$Sample_ID %in% colnames(Rcountdata))
dds <- DESeqDataSetFromMatrix(countData = Rcountdata,
                              colData = Rcovariates, 
                              design = ~ Dex)

#Filter out genes that are not expressed in this brain region
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]


# 3. Normalization, box-plots and MA-plots ----
#Sequencing depth normalization
dds <- estimateSizeFactors(dds)

# Extract counts and size factors
norm_counts <- counts(dds, normalized = TRUE)
raw_counts <- counts(dds, normalized = FALSE)


# Box-plots before and after the normalization
png(filename = paste0(prefix_plots, "readCounts_nonnormalized.png"))
bp_before <- boxplot(log2(counts(dds)+1), notch=TRUE,
                     main = "Non-normalized read counts",
                     ylab="log2(read counts)", cex=.6)
dev.off()

png(filename = paste0(prefix_plots, "readCounts_normalized.png"))
bp_after <- boxplot(log2(counts(dds, normalize = TRUE)+1), notch=TRUE,
                    main = "Size-factor-normalized read counts",
                    ylab="log2(read counts)", cex=.6)
dev.off()


# 4. Unsupervised Clustering Analysis and Batch Correction ----
# 4.1 Transformation ----
# 4.1.1 VST
# Variance Stabilizing Transformation (VST) - log transformation moderates the variance accross the mean
vsd <- vst(dds, blind = TRUE)

# 4.1.2 compare log2+1 and vst by plotting first sample against second for each method
df <- bind_rows(
  as.data.frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)


# 4.2 Sample correlations -------------
vsd_cor_val <- cor(assay(vsd))
# Observe the sample distances
vsd_cor_val_group <- dplyr::select(Rcovariates, Dex)
SampleOrder <- order(vsd_cor_val_group$Dex)
pheatmap(vsd_cor_val[,SampleOrder], 
         annotation_col = vsd_cor_val_group)
pheatmap(vsd_cor_val[,SampleOrder], 
         annotation_col = vsd_cor_val_group,
         filename = paste0(prefix_plots, "samplecorrelations.png"))
dev.off()


# 5. PCA analysis and plots ----------------
# calculate PCs and add to summary table
pc <- pca(vsd@assays@data[[1]], removeVar = 0.2)
colData(dds) <- DataFrame(cbind(colData(dds), pc$rotated[,c(1:10)] ))
head(colData(dds))

# PCA plots
percentVar <- round(pc$variance[1:10], digits = 1)
# Dex
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = Dex)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_Dex.png"))
# possible noise factors
# mouse weight
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = `mouse_weight.g.`)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_MouseWeight.png"))
# injection volume
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = `Injection_volume`)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_InjectionVolume.png"))
# cryostat machine
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = cryostat)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_Cryostat.png"))
# researcher
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = `Researcher`)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_Researcher.png"))
# date of punching
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = `date_of_punching`)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_DateOfPunching.png"))
# plate, lane, row
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = Plate)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_Plate.png"))
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = Lane, shape = Plate)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_LanePlate.png"))
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = Row, shape = Plate)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_RowPlate.png"))
# RIN and RIN2
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = RIN)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_RIN.png"))
ggplot(as.data.frame(colData(dds)), aes(x = PC1, y = PC2, color = RIN2)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")
ggsave(paste0(prefix_plots, "pca_RIN2.png"))

#pairs plot
#Color Dex
png(filename = paste0(prefix_plots, "pcaPairsplot_Dex.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$Dex,
      cex=2,pch='*'))
dev.off()
#Color date
png(filename = paste0(prefix_plots, "pcaPairsplot_DateOfPunching.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$date_of_punching,
      cex=2,pch='*'))
dev.off()
#Color Researcher
png(filename = paste0(prefix_plots, "pcaPairsplot_Researcher.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$Researcher,
      cex=2,pch='*'))
dev.off()
#Color Cryostat
png(filename = paste0(prefix_plots, "pcaPairsplot_Cryostat.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$cryostat,
      cex=2,pch='*'))
dev.off()
#Color Plate
png(filename = paste0(prefix_plots, "pcaPairsplot_Plate.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$Plate,
      cex=2,pch='*'))
dev.off()
#Color Lane
png(filename = paste0(prefix_plots, "pcaPairsplot_Lane.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$Lane,
      cex=2,pch='*'))
dev.off()
#Color Row
png(filename = paste0(prefix_plots, "pcaPairsplot_Row.png"), width = 900, height = 900)
print(splom(as.data.frame(pc$rotated[,1:10]),
      col=colData(dds)$Row,
      cex=2,pch='*'))
dev.off()


#MDS plot (not necessary if we have the complete matrix, helpful if only distances are known)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Dex, label = Sample_ID)) +
  geom_point(size = 3) + 
  geom_text_repel() + 
  coord_fixed() + ggtitle("MDS with VST data")
ggsave(paste0(prefix_plots, "mds_Dex.png"))


# 6. Outlier detection on dex and baseline samples separately ---------------
dds <- outlier_removal(dds)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)


# 7. Batch Identification with SVA (SVA,CCA) ----------------
batch <- batch_sva(dds)
cov_pc <- batch$cov_pc
n.sv <- batch$n.sv


# 8. Corrected model ---------------
# 8.1 Add SVs to model design -----
ddssva <- dds
colData(ddssva) <- cbind(colData(ddssva), cov_pc[,paste0("SV", seq(1:n.sv))])
design(ddssva) <- as.formula(paste0("~ Dex +", paste(paste0("SV", seq(1:n.sv)), collapse = "+")))


# 8.2 Variance Partition ----------
form <- as.formula(paste0("~ ", paste0( grep(colnames(colData(ddssva)), pattern="SV", value=T), collapse=" + " ), 
                          "+ (1|Dex)", collapse=""))
form

# run variance partition
varPart <- fitExtractVarPartModel(assay(vsd), form, as.data.frame(colData(ddssva)))

# plot variance partition (Violin plot of variance fraction for each gene and each variable)
png(paste0(prefix_plots, "variancePartition.png"))
plotVarPart(varPart)
dev.off()

# plot percentage of variance explained by each variable
print("% variance explained by covariates:")
o <- colSums(varPart)*100/sum(varPart)
o
png(paste0(prefix_plots, "varPart_percentage.png"))
barplot(o, las=2)
dev.off()



# 9. Differential Expression Analysis ------------
# 9.1 Running DE analysis based on the Negative Binomial (aka Gamma-Poisson) distribution
# DESeq method performs three steps:
# a. estimation of size factors (controlling for differences in the sequencing depth)
# b. estimation of dispersion values for each gene
# c. fitting a generalized linear model
# ddssva$Dex <- relevel(ddssva$Dex, ref = "1")
ddssva <- DESeq(ddssva)
plotDispEsts(ddssva)

# 9.2 Building the results table 
# extract estimated log2 fold changes and p values 
# without contrast: default is last variable in design
res <- results(ddssva, contrast=c("Dex","1","0"))
# positive logFC: high in dex, low in control using this contrast
# negative logFC: high in control, low in dex treated samples
res
summary(res)
# information on the meaning of the columns in res
mcols(res, use.names = TRUE)
# The first column, baseMean, is a just the average of the normalized count values, 
# divided by the size factors, taken over all samples in the DESeqDataSet. 
# The remaining four columns refer to a specific contrast, namely the comparison of 
# the trt level over the untrt level for the factor variable dex. We will find out 
# below how to obtain other contrasts.
# The column log2FoldChange is the effect size estimate. It tells us how much the 
# gene’s expression seems to have changed due to treatment with dexamethasone in 
# comparison to untreated samples. This value is reported on a logarithmic scale to 
# base 2: for example, a log2 fold change of 1.5 means that the gene’s expression is 
# increased by a multiplicative factor of 21.5≈2.82.
# Of course, this estimate has an uncertainty associated with it, which is available 
# in the column lfcSE, the standard error estimate for the log2 fold change estimate. 
# We can also express the uncertainty of a particular effect size estimate as the 
# result of a statistical test. The purpose of a test for differential expression 
# is to test whether the data provides sufficient evidence to conclude that this 
# value is really different from zero. DESeq2 performs for each gene a hypothesis 
# test to see whether evidence is sufficient to decide against the null hypothesis 
# that there is zero effect of the treatment on the gene and that the observed 
# difference between treatment and control was merely caused by experimental variability 
# (i.e., the type of variability that you can expect between different samples in the 
#   same treatment group). As usual in statistics, the result of this test is reported 
# as a p value, and it is found in the column pvalue. Remember that a p value indicates 
# the probability that a fold change as strong as the observed one, or even stronger, 
# would be seen under the situation described by the null hypothesis.

# lower the FDR in result table
res.05 <- results(ddssva, contrast=c("Dex","1","0"), alpha = 0.05)
table(res.05$padj < 0.05)

# raise the logFC threshold from 0
resLFC1 <- results(ddssva, contrast=c("Dex","1","0"), lfcThreshold=1)
table(resLFC1$padj < 0.1)

# #plot gene with smallest pvalue
# topGene <- rownames(res)[which.min(res$padj)]
# plotCounts(ddssva, gene = topGene, intgroup=c("Dex"))
# library("ggbeeswarm")
# geneCounts <- plotCounts(ddssva, gene = topGene, intgroup = c("Dex","Mouse_ID"),
#                          returnData = TRUE)
# ggplot(geneCounts, aes(x = Dex, y = count, color = Mouse_ID)) +
#   scale_y_log10() +  geom_beeswarm(cex = 3)


# 9.3 Shrink log2 fold changes
# In statistics, shrinkage is the reduction in the effects of sampling variation. 
# In regression analysis, a fitted relationship appears to perform less well on a 
# new data set than on the data set used for fitting.[1] In particular the value of 
# the coefficient of determination 'shrinks'. This idea is complementary to overfitting 
# and, separately, to the standard adjustment made in the coefficient of determination 
# to compensate for the subjunctive effects of further sampling, like controlling for 
# the potential of new explanatory terms improving the model by chance: that is, the 
# adjustment formula itself provides "shrinkage." But the adjustment formula yields 
# an artificial shrinkage.
# A shrinkage estimator is an estimator that, either explicitly or implicitly, incorporates 
# the effects of shrinkage. In loose terms this means that a naive or raw estimate is 
# improved by combining it with other information. The term relates to the notion that 
# the improved estimate is made closer to the value supplied by the 'other information' 
# than the raw estimate. In this sense, shrinkage is used to regularize ill-posed 
# inference problems. (Source: Wikipedia)
resultsNames(ddssva)
res <- lfcShrink(ddssva, coef="Dex_1_vs_0", type = "apeglm")
res <- res[order(res$padj),]
head(res)
res_print <- as.data.table(res, keep.rownames = TRUE) %>%
  dplyr::rename("Ensembl_ID" = "rn")
res_print$Gene_Symbol <- mapIds(org.Mm.eg.db, keys = res_print$Ensembl_ID, 
                          keytype = "ENSEMBL", column="SYMBOL")
head(res_print)
write.table(res_print, file=paste0(prefix_tables, "Dex_1_vs_0_lfcShrink.txt"),
            sep="\t", quote = F, row.names = FALSE)

# 9.4 Plots
# MA plot (mean of normalized counts vs LFC)
png(paste0(prefix_plots, "MAplot_lfcShrink.png"))
DESeq2::plotMA(res,ylim = c(-5, 5))
dev.off()
# MA plot without shrinkage
res.noshr <- results(ddssva, name="Dex_1_vs_0")
DESeq2::plotMA(res.noshr, ylim = c(-5, 5))

# volcano plot
png(paste0(prefix_plots, "volcanoplot_lfcShrink.png"),
    width = 600,
    height = 600)
print(EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "pvalue",
                title = paste0("DESeq2 results: ", region),
                subtitle = "Differential expression"))
dev.off()

volcano_data <- as.data.frame(res@listData)
volcano_data$sig <- as.factor(volcano_data$padj <= 0.1)
ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=sig)) +
  geom_point() + 
  scale_color_manual(values = c("#FF9900", "lightgrey"),
                     breaks = c("TRUE", "FALSE"),
                     labels = c("padj <= 0.1", "padj > 0.1"),
                     name = "") +
  xlab("log2(FoldChange)") +
  ylab("-log10(padj)") +
  theme_bw() +
  theme(text = element_text(size= 20)) +
        #legend.position = "bottomright") +
  labs(title = element_text("Differential expression analysis \nwith DESeq2"))
ggsave(filename = paste0(prefix_plots, "volcanoplot_lfcShrink.svg"),
       width = 8, height = 6)

# # plot gene with smallest p value
# DESeq2::plotMA(res,ylim = c(-5, 5))
# topGene <- rownames(res)[which.min(res$padj)]
# with(res[topGene, ], {
#   points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
#   text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
# })
# 
# # plot histogram of pvalue distribution
# hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
#      col = "grey50", border = "white")


# # plot gene expression of most variable genes as heatmap
# library("genefilter")
# topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
# mat  <- assay(vsd)[ topVarGenes, ]
# mat  <- mat - rowMeans(mat)
# anno <- as.data.frame(colData(vsd)[, "Dex"])
# rownames(anno) <- colnames(mat)
# pheatmap(mat, annotation_col = anno)




# 10. Differential Expression Analysis with lfcThreshold (analog to treat function in limma) ------------
# 10.1 Running DE analysis based on the Negative Binomial (aka Gamma-Poisson) distribution
#ddssva <- DESeq(ddssva)
#plotDispEsts(ddssva)

# 10.2 Building the results table 
# extract estimated log2 fold changes and p values 
res <- results(ddssva, contrast=c("Dex","1","0"), lfcThreshold = 0.5)
# downregulation (negative logFC): high in dex 0, low in dex 1 using this contrast
res <- res[order(res$padj),]
res
summary(res)
# information on the meaning of the columns in res
mcols(res, use.names = TRUE)


# 10.3 Shrink log2 fold changes
# svalues: The adjusted p-values and s-values are similar but with a different 
# definition of error. One focuses on falsely rejecting what are truly null genes, 
# and the other on getting the sign of the LFC wrong.
res_shrink <- lfcShrink(ddssva, coef = "Dex_1_vs_0", type = "apeglm", lfcThreshold = 0.5)
res_shrink <- res_shrink[order(res_shrink$svalue),]
head(res_shrink)
res_shrink_print <- as.data.table(res_shrink, keep.rownames = TRUE) %>%
  dplyr::rename("Ensembl_ID" = "rn")
res_shrink_print$Gene_Symbol <- mapIds(org.Mm.eg.db, keys = res_shrink_print$Ensembl_ID, 
                                keytype = "ENSEMBL", column="SYMBOL")
head(res_shrink_print)
write.table(res_shrink_print, file=paste0(prefix_tables, "Dex_1_vs_0_lfcShrink_lfc0.5.txt"),
            sep="\t", quote = F, row.names = FALSE)

# 10.4 Plots
# MA plot (mean of normalized counts vs LFC)
png(paste0(prefix_plots, "MAplot_lfcShrink_lfc0.5.png"))
DESeq2::plotMA(res_shrink,ylim = c(-5, 5))
dev.off()
# MA plot without shrinkage
DESeq2::plotMA(res, ylim = c(-5, 5))

# volcano plot
png(paste0(prefix_plots, "volcanoplot_lfcShrink_lfc0.5.png"),
    width = 600,
    height = 600)
print(EnhancedVolcano(res_shrink,
                lab = rownames(res_shrink),
                x = "log2FoldChange",
                y = "svalue",
                title = paste0("DESeq2 results: ", region),
                subtitle = "Differential expression"))
dev.off()
# volcano plot without shrinkage
EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "pvalue",
                title = paste0("DESeq2 results: ", region),
                subtitle = "Differential expression")


# 11. Print data in correct format for kimono ------------
print_kimono(vsd, colData(ddssva))
