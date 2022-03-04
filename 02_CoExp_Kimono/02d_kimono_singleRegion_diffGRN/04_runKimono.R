#!/usr/bin/env Rscript

###############################################

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
# if (length(args)!=5) {
#   stop("Please supply expression.csv, phenotypes.csv, prior_expr_bio.csv and output file", call.=FALSE)
# }
reg <- "PVN"
d <- 0

##############################################
# manual(if not snakemake): data files
args[[1]]=paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/expression_mm_",reg,"_dex",d,".csv")
args[[2]]=paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/phenotypes_mm_",reg,"_dex",d,".csv")
args[[3]]="/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/prior_expr_bio_mm_regDex.csv"
args[[4]]="/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/data/kimono_input/prior_expr_funcoup_mm.csv"
args[[5]]=paste0("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/tables/coExpression_kimono/04_singleRegion_",reg,"_dex",d,"_funcoup.csv")

###############################################
### 1 libraries & code
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify", "stringr",
               "magrittr", "foreach", "doParallel", "dplyr", "furrr", "purrr")
lapply(libraries, require, character.only = TRUE)

source("/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/scripts/02_CoExp_Kimono/kimono_stability/infer_sgl_model.R")
source('/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/scripts/02_CoExp_Kimono/kimono_stability/utility_functions.R')
source('/Users/nathalie_gerstner/Documents/ownCloud/DexStim_RNAseq_Mouse/scripts/02_CoExp_Kimono/kimono_stability/kimono.R')

print("libraries & moni main done")



###############################################
### 2 Data
###############################################

mrna <- as.data.frame(as.matrix(fread(args[1]),rownames=1))
# mrna <- mrna[1:5,]


###############################################

biological <- as.data.frame(as.matrix(fread(args[2]),rownames=1))

# (make sure the rows are in the same order)
idorder <- as.character(rownames(biological))
mrna <- mrna[match(idorder, rownames(mrna)),]


print("read input")
###############################################
###############################################

# mapping

prior_mrna_bio <- fread(args[[3]]);prior_mrna_bio[1:5,]
prior_mrna <- fread(args[[4]]); prior_mrna[1:5,]

print("read mapping")




###############################################
# 3 Assemble into lists
###############################################

#set input parameters
input_list <- list(
  as.data.table(mrna),
  # as.data.table(mrna),
  as.data.table(biological)
)
names(input_list) <- c('mrna',
                       'biological')

#########################
mapping_list <- list(
  as.data.table(prior_mrna_bio),
  as.data.table(prior_mrna)
)
#########################
metainfo <-data.frame('ID'   = c('mrna_bio', 'prior_mrna'),
                      'main_to'   =  c(2,1)
)
print("data created")


###############################################
# 4 Run MONI
###############################################

# parallel
options(future.globals.maxSize = 30000 * 1024^2) # more memory for each thread
plan(multisession, workers = 12) # 12 parallel



# start
start_time <- Sys.time();start_time

node_list <- colnames(input_list$mrna)[1:10]
results <- future_map(node_list, run_kimono_para, stab_sel = TRUE, .progress = TRUE, .options = future_options(seed = TRUE))

# falls nicht parallelisiert:
# results <-  kimono(input_list, mapping_list, metainfo, main_layer = 1, min_features = 5, sel_iterations = 10, core = 2)


Sys.time()- start_time
results <- do.call(rbind, results) # make one big data table

###############################################


print("kimono done"); end_time <- Sys.time();end_time - start_time



#end

fwrite(results,file=args[5],
       row.names = FALSE,quote = FALSE)
print("table written")


