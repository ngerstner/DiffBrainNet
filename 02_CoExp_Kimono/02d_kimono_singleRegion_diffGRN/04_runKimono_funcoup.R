#!/usr/bin/env Rscript

###############################################

args = commandArgs(trailingOnly=TRUE)

# read region, dex status and startnode
reg <- args[[1]]
#reg <- "dCA1"
d <- args[[2]]
#d <- 0
startnode <- as.numeric(args[[3]])

##############################################
# manual(if not snakemake): data files
basepath <- "/binder/mgp/workspace/DexStim_RNAseq_Mouse/"
input_expr <- paste0(basepath,"data/kimono_input/expression_mm_",reg,"_dex",d,".csv")
input_pheno <- paste0(basepath,"data/kimono_input/phenotypes_mm_",reg,"_dex",d,".csv")
input_prior_bio <- paste0(basepath,"data/kimono_input/prior_expr_bio_mm_regDex.csv")
input_prior <- paste0(basepath,"data/kimono_input/prior_expr_funcoup_mm.csv")
output <- paste0(basepath,"tables/coExpression_kimono/04_singleRegion_",reg,"_dex",d,"_funcoup_parallel_",startnode,".csv")

###############################################
### 1 libraries & code
###############################################

#.libPaths( c( .libPaths(), "/binder/mgp/workspace/DexStim_RNAseq_Mouse/Rpackages/") )
#Sys.setenv(R_LIBS = paste("/binder/mgp/workspace/DexStim_RNAseq_Mouse/Rpackages/", Sys.getenv("R_LIBS"), sep=.Platform$path.sep))
libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify", "stringr",
               "magrittr", "foreach", "doParallel", "dplyr", "furrr", "purrr")
lapply(libraries, require, character.only = TRUE)

source(paste0(basepath,"scripts/02_CoExp_Kimono/kimono_stability/infer_sgl_model.R"))
source(paste0(basepath,'scripts/02_CoExp_Kimono/kimono_stability/utility_functions.R'))
source(paste0(basepath,'scripts/02_CoExp_Kimono/kimono_stability/kimono.R'))

print("libraries & moni main done")



###############################################
### 2 Data
###############################################

mrna <- as.data.frame(as.matrix(fread(input_expr),rownames=1))
# mrna <- mrna[1:5,]


###############################################

biological <- as.data.frame(as.matrix(fread(input_pheno),rownames=1))

# (make sure the rows are in the same order)
idorder <- as.character(rownames(biological))
mrna <- mrna[match(idorder, rownames(mrna)),]


print("read input")
###############################################
###############################################

# mapping

prior_mrna_bio <- fread(input_prior_bio);prior_mrna_bio[1:5,]
prior_mrna <- fread(input_prior); prior_mrna[1:5,]

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
plan(multisession, workers = 12)

# start
start_time <- Sys.time();start_time


endnode <- min(startnode+1000, length(colnames(input_list$mrna)))
node_list <- colnames(input_list$mrna)[startnode:endnode]
results <- future_map(node_list, run_kimono_para, stab_sel = TRUE, niterations = 20, .options = furrr_options(seed = TRUE)) #add one more layer of parallelization in kimono.R (seeds)


Sys.time()- start_time
results <- do.call(rbind, results) # make one big data table

###############################################


print("kimono done"); end_time <- Sys.time();end_time - start_time



#end

fwrite(results,file=output,
       row.names = FALSE,quote = FALSE)
print("table written")


