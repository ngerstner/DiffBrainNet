##################################################
## Project: DexStim Mouse Brain
## Date: 22.07.2020
## Author: Nathalie
###################################################
# RNA-seq Analysis. Functions that are needed more often

# Outlier detection on dex and baseline samples separately
outlier_removal <- function(dds) {
  # Dex samples
  outlier <- TRUE
  iter <- 1
  while (outlier) {
    # Run PCA
    dds_dex <- dds[, colData(dds)$Dex == 1]
    dds_dex <- estimateSizeFactors(dds_dex)
    vsd_dex <- vst(dds_dex, blind = TRUE)
    pc_dex <-
      as.data.frame(pca(vsd_dex@assays@data[[1]], removeVar = 0.2)$rotated)
    pc_dex$sample <- rownames(pc_dex)
    splom(as.data.frame(pc_dex[, 1:10]),
          cex = 2, pch = '*')
    
    # Label any sample which is more than 2.5SD away from the mean in PC1 as outlier
    outlier_dex <-
      which(abs(pc_dex[, 1] - mean(pc_dex[, 1])) > (2.5 * sd(pc_dex[, 1])))
    pc_dex$outlier <- FALSE
    pc_dex$outlier[outlier_dex] <- TRUE
    
    ggplot(pc_dex, aes(
      x = PC1,
      y = PC2,
      color = outlier,
      label = sample
    )) +
      geom_point(size = 3) +
      geom_text_repel() +
      xlab("PC1") +
      ylab("PC2") +
      coord_fixed() +
      ggtitle("PCA: Dex Samples")
    ggsave(paste0(prefix_plots, "outlier_Dex_iter", iter, ".png"))
    
    # Subset deseq object / remove outliers
    keep_names <-
      setdiff(colData(dds)$Sample_ID, pc_dex$sample[outlier_dex])
    dds <- dds[, keep_names]
    iter <- iter + 1
    if (length(outlier_dex) == 0) {
      outlier <- FALSE
    }
  }
  
  # Baseline samples
  outlier <- TRUE
  iter <- 1
  while (outlier) {
    # Run PCA
    dds_base <- dds[, colData(dds)$Dex == 0]
    dds_base <- estimateSizeFactors(dds_base)
    vsd_base <- vst(dds_base, blind = TRUE)
    pc_base <-
      as.data.frame(pca(vsd_base@assays@data[[1]], removeVar = 0.2)$rotated)
    pc_base$sample <- rownames(pc_base)
    splom(as.data.frame(pc_base[, 1:10]),
          cex = 2, pch = '*')
    
    # Label any sample which is more than 2.5SD away from the mean in PC1 as outlier
    outlier_base <-
      which(abs(pc_base[, 1] - mean(pc_base[, 1])) > (2.5 * sd(pc_base[, 1])))
    pc_base$outlier <- FALSE
    pc_base$outlier[outlier_base] <- TRUE
    
    ggplot(pc_base, aes(
      x = PC1,
      y = PC2,
      color = outlier,
      label = sample
    )) +
      geom_point(size = 3) +
      geom_text_repel() +
      xlab("PC1") +
      ylab("PC2") +
      coord_fixed() +
      ggtitle("PCA: Baseline Samples")
    ggsave(paste0(prefix_plots, "outlier_Baseline_iter", iter, ".png"))
    
    # Subset deseq object / remove outliers
    keep_names <-
      setdiff(colData(dds)$Sample_ID, pc_base$sample[outlier_base])
    dds <- dds[, keep_names]
    iter <- iter + 1 
    if (length(outlier_base) == 0) {
      outlier <- FALSE
    }
  }
  
  # Rerun normalization etc on combined dataset
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE)
  pc <-
    as.data.frame(pca(vsd@assays@data[[1]], removeVar = 0.2)$rotated)
  pc$Dex <- colData(dds)$Dex
  
  ggplot(pc, aes(x = PC1, y = PC2, color = Dex)) +
    geom_point(size = 3) +
    xlab("PC1") +
    ylab("PC2") +
    coord_fixed() +
    ggtitle("PCA: Outlier Deleted")
  ggsave(paste0(prefix_plots, "outlierDeleted.png"))
  
  png(
    filename = paste0(prefix_plots, "pairsplot_outlierDeleted.png"),
    width = 900,
    height = 900
  )
  print(splom(
    as.data.frame(pc[, 1:10]),
    col = colData(dds)$Dex,
    cex = 2,
    pch = '*'
  ))
  dev.off()
  
  # add new PCs to colData
  colData(dds)[, c(26:35)] <- pc[, c(1:10)]
  
  return(dds)
}



# Batch identification with SVA
batch_sva <- function(dds){
  #possible batch effects
  #Explanations: Sample_ID is not a batch, Animal is same as Mouse_ID, mouse_weight.g. is linearly correlated 
  #with Injection_volume, also colinearity between Researcher, date_of_punching and cryostat,
  #sample wells and indices are covered by plate. lane and row ...
  # --> maybe remove conc.RNA.ng.ml and vol.for.25.ng.input --> colinearity to Dex? --> should not be removed
  pos_batch <- c("Dex", "Injection_volume", "date_of_punching",
                 "Plate", "Lane", "Row", "RIN", "RIN2")
  cov_pc <- colData(dds)[,c(pos_batch,paste0("PC", seq(1:10)))]
  
  
  # 1. Surrogate Variable Analysis (SVA) ----
  mod <- model.matrix(~ Dex, colData(dds))
  mod0 <- model.matrix(~ 1, colData(dds))
  # Calculate SVs (on normalized counts --> according to documentation)
  # Comment: don't use num.sv function, apparently it's only for microarray data
  norm <- counts(dds, normalized = TRUE)
  svobj <- svaseq(norm, mod, mod0)
  n.sv <- svobj$n.sv #Number of significant surrogate variables is
  # Add significant SVs to covariates
  coln <- colnames(cov_pc)
  cov_pc <- cbind(cov_pc,svobj$sv[,1:n.sv])
  colnames(cov_pc) <- c(coln, paste0("SV", seq(1:n.sv)))
  
  
  # 2. Canonical Correlation Analysis (CCA) ----
  form <- as.formula(paste0("~ Dex +
  Injection_volume +
  date_of_punching +
  Plate + Row + Lane +
  RIN + RIN2 +
  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+",
  paste(paste0("SV", seq(1:n.sv)), collapse = "+")))
  
  # Calculate the correlation coefficients
  C <- canCorPairs(form, cov_pc)
  # Plot the results using Canonical correlation
  png(filename = paste0(prefix_plots, "cca_sorted.png"), width = 800, height = 800)
  plotCorrMatrix(C)
  dev.off()
  png(filename = paste0(prefix_plots, "cca_unsorted.png"), width = 800, height = 800)
  plotCorrMatrix(C, sort = FALSE)
  dev.off()
  
  
  # 3. P-values for correlations
  pval_corr <- matrix(1, nrow = 10 + n.sv, ncol = ncol(cov_pc) - (10 + n.sv))
  rownames(pval_corr) <-
    c(paste0("PC", seq(1:10)), paste0("SV", seq(1:n.sv)))
  colnames(pval_corr) <- pos_batch
  # Calc pvalue for factor covariates using ANOVA
  for (cov in names(Filter(is.factor, cov_pc))){
    for (v in c(paste0("PC", seq(1:10)), paste0("SV", seq(1:n.sv)))) {
      f <- paste0(v, "~", cov)
      p <- summary(aov(as.formula(f), data = cov_pc))[[1]]$`Pr(>F)`[[1]]
      pval_corr[v, cov] <- p
    }
  }
  # Calc pvalues for numeric covariates using linear model
  for (cov in names(Filter(is.numeric, cov_pc[,1:(ncol(cov_pc)-(10+n.sv))]))){
    for (v in c(paste0("PC", seq(1:10)), paste0("SV", seq(1:n.sv)))) {
      f <- paste0(v, "~", cov)
      p <- summary(lm(as.formula(f), data = cov_pc))$coefficients[2,4]
      pval_corr[v, cov] <- p
    }
  }
  
  # Plot p-values
  pval_corr
  pheatmap(pval_corr, cluster_rows = FALSE)
  pheatmap(pval_corr)
  pheatmap(pval_corr, cluster_rows = FALSE,
           filename = paste0(prefix_plots, "batchevaluation_pval.png"))
  dev.off()
  
  return(list("cov_pc" = cov_pc, "n.sv" = n.sv))
}


# Print data in correct format for kimono (at least for this region)
print_kimono <- function(vsd, cov_data){
  
  # Write vst transformed data to file
  write.table(assay(vsd), file=paste0(prefix_tables, "expression_vsd.txt"),sep="\t", quote = F)
  
  # Write biological variables and SVs to file
  biol <- cov_data[,c("Dex", "Region", grep(colnames(colData(ddssva)), pattern="SV", value=T))]
  write.table(biol, file=paste0(prefix_tables, "bio_variables.txt"),sep="\t", quote = F)
}