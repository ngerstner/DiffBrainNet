##################################################
## Project: DexStim Mouse Brain
## Date: 22.07.2020
## Author: Nathalie
###################################################
# RNA-seq Analysis. Part 1. Setting up the environment

# This R script "setup" is for:
# 1. loading the nessesary packages
# 2. loading and cleaning data (metha and readcounts)
# 3. pre-filtering
# 4. perfom simple descriptive statistics

# set working directory to source file location
setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

# 1. Packages ----

## Library Loadings
library(DESeq2) # DE analysis
library(limma) # DE analysis (and MAplots for DESeq2)
library(sva) # Surrogate variable analysis (SVA)
library(Biobase)
library(variancePartition) # canonical correlation analysis (CCA)

# 2. Loading and cleaning the data ----

# Data pathways
db_metha <- "data/Covariates_RNAseq_Mouse_brain.csv"
db_data <- "data/06_featureCounts/20190313_Anthi_Mouse_Brain_Dex.fC"

# Loading
covariates <- read.csv(db_metha)
countdata <- read.csv(file = db_data, sep = "\t", header = TRUE, row.names = 1) # removed very first line before actual header
countdata <- countdata[,6:ncol(countdata)] #Remove first six columns (geneid, chr, start, end, strand, length)
ncol(countdata) #[1] 297
nrow(countdata) #[1] 27998

# Adjust the sample names
colnames(countdata) <- gsub("AK_S[0-9]*.cutadapt.extract.Aligned.sortedByCoord.out.dedup.bam", "", colnames(countdata))
colnames(countdata) <- gsub("X05_umi_dedup.mpg_L[0-9]*_", "", colnames(countdata))
colnames(countdata) <- sub("[.]", "_", colnames(countdata))
colnames(countdata) <- sub("[.]", "", colnames(countdata))

# Adjust the data in files (delete unknown or not-matching)
all_col_names <- colnames(countdata)
cov_row_names <- as.vector(covariates$Sample_ID)

outersect <- function(x, y) { #function to identify the non-shared items in two vectors
  sort(c(setdiff(x, y), setdiff(y, x)))}
to_delete <- outersect(all_col_names, cov_row_names) 
print(to_delete)

covariates <- covariates[!covariates$Sample_ID %in% to_delete,]
countdata <- countdata[,!colnames(countdata) %in% to_delete]

# 3. Pre-filtering ----

# Convert to matrix
countdata <- as.matrix(countdata)

# 3.1 Gene-based filtering ---- 
# Keep genes quantified in at least one full treatment group
region <- as.character(unique(covariates$Region))
dex <- unique(covariates$Dex)
filtered_data <- matrix(nrow = nrow(countdata), ncol = ncol(countdata))

# Go through all combinations of brain region and dex treatment
for (r in region){
  for (d in dex){
    sample_ids <- covariates$Sample_ID[covariates$Region == r & covariates$Dex == d]
    # Subset counts to samples of this dex/region group
    sub_counts <- countdata[,colnames(countdata)%in%sample_ids]
    # Copy gene expression data if a gene is expressed in all samples of the dex/region group
    for (k in 1:nrow(countdata)) { 
      if (length(which(sub_counts[k,]==0)) == 0){
        filtered_data[k,] <- countdata[k,]
      }
    }
  }
}
colnames(filtered_data) <- colnames(countdata)
rownames(filtered_data) <- rownames(countdata)
filtered_data <- na.omit(filtered_data)
# 12976 genes/rows and 295 samples/columns 


# 3.2 Sample-based filtering ---- 
# Check 0 read counts (assess if there are samples with way too low counts --> indicating errors)
# Should be more than 60%

# 90% of the data present is roughly 11700 => 1300 may be zeros
#With this cut-off there are 291 samples left (4 are deleted)
# 80% of the data present is roughly 10400 => 2600 may be zeros
#With this cut-off there are 295 samples left (none are deleted)
# 70% of the data present is roughly 9100 => 3900 may be zeros
#With this cut-off there are 295 samples left (none are deleted)

delete_ids <- NULL

for (i in 1:ncol(filtered_data)){
  if (length(which(filtered_data[,i]==0)) > 3900){
    delete_ids <- c(delete_ids, colnames(filtered_data)[i])
  }
}

filtered_data <- filtered_data[,!colnames(filtered_data) %in% delete_ids]


# 4. Simple descriptive statistics (Iuliia) ----

# Genes per brain region/condition table (data$Region, data$Dex)
## To report how many genes are left after pre-filtering
genes_per_reg_cond <- matrix(nrow = length(region), ncol = length(dex))
rownames(genes_per_reg_cond) <- region
colnames(genes_per_reg_cond) <- dex
for (r in region){
  for (d in dex){
    sample_ids <- covariates$Sample_ID[covariates$Region == r & covariates$Dex == d]
    sub_counts <- filtered_data[,colnames(filtered_data) %in% sample_ids]
    sub_counts <- rbind(sub_counts, NA)
    #keep the present genes number and calculate the median value of those
    for (k in 1:ncol(sub_counts)) { 
      sub_counts[nrow(sub_counts), k] <- length(which(sub_counts[,k] > 0))
    }
    genes_per_reg_cond[r, toString(d)] <- median(as.numeric(sub_counts[nrow(sub_counts),]))
  }
}

## Adjust the tables again
#####!!!!!Check why countdata delets here
all_col_names <- colnames(countdata)
sample_filter_col_names <- colnames(filtered_data)
to_delete <- outersect(all_col_names,sample_filter_col_names)
covariates <- as.data.frame(covariates[!covariates$Sample_ID %in% to_delete,]) #data.frame for DESeq2
countdata <- (countdata[,!colnames(countdata) %in% to_delete])

data_row_names <- rownames(countdata)
sample_filter_row_names <- rownames(filtered_data)
to_delete <- outersect(data_row_names, sample_filter_row_names)
countdata <- as.matrix(countdata[!rownames(countdata) %in% to_delete,])
backupcov<-covariates
#ncol(countdata) #295
#nrow(countdata) #12976

#Build a histogram of mouse weight
hist(as.numeric(covariates$mouse_weight.g.[covariates$Dex==0]))
hist(as.numeric(covariates$mouse_weight.g.[covariates$Dex==1]))
shapiro.test(as.numeric(covariates$mouse_weight.g.[covariates$Dex==0]))     
shapiro.test(as.numeric(covariates$mouse_weight.g.[covariates$Dex==1]))     
wilcox.test(as.numeric(covariates$mouse_weight.g.)~covariates$Dex)
ggplot(covariates,aes(x=mouse_weight.g., color=Dex, fill=Dex))+
  geom_histogram(alpha=0.6,  binwidth = 0.5)


# Assign Conditions for the Data Analysis ----
str(covariates)
covariates$Dex = as.factor(covariates$Dex)
covariates$Region = as.factor(covariates$Region)
covariates$Region <- relevel(covariates$Region, "CER")
covariates$Animal = as.factor(covariates$Animal)
covariates$Injection_volume <- as.factor(covariates$Injection_volume)
covariates$Researcher = as.factor(covariates$Researcher)
covariates$date_of_punching = as.factor(covariates$date_of_punching)
covariates$Plate = as.factor(covariates$Plate)
covariates$Lane = as.factor(covariates$Lane)
covariates$Row = as.factor(covariates$Row)
covariates$cryostat = as.factor(covariates$cryostat)
covariates$Sample_Well = as.factor(covariates$Sample_Well)
covariates$index = as.factor(covariates$index)
covariates$index2 = as.factor(covariates$index2)
covariates$I5_Index_ID = as.factor(covariates$I5_Index_ID)
covariates$I7_Index_ID = as.factor(covariates$I7_Index_ID)

# sort covariates table to have same order as countdata
rownames(covariates) <- covariates$Sample_ID
idx<-match(colnames(countdata),rownames(covariates))
covariates<-covariates[idx,]

