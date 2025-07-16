#!/usr/bin/env Rscript
### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# getwd()
gc()
memory.limit(size=Inf)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, data.table, foreach, doParallel, preprocessCore) # parallel

parser <- OptionParser()
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c( "--gene_train"), default="./data/atlas_gene.tsv",
              help = "Set file path for Gene Expression data used for training [default \"%default\"]. In tab-delimited format.",
              metavar="training gene expression"),
  make_option(c( "--gene_test"), default="./data/atlas_gene.tsv",
              help = "Set file path for Gene Expression data used for testing [default \"%default\"]. In tab-delimited format.",
              metavar="training gene expression"),  
  make_option(c("--proj_name"), default="new_proj",
              help = "Define a project name. Distance metrics will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--dist"), default="pcc",
              help = "Method to use to calculate nearest training cell-types. 'pcc' for Pearson Correlation, 'scc' for Spearman Correlation, and 'euc' for Euclidean distance)",
              metavar="dist"),
  make_option(c("--task_id"), default="1", help = "Task identifier (user-specified) for this experiment", metavar = "task id")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(opt$verbose){
  print(opt)
}

wd <- opt$wd
setwd(wd)

dir.create(sprintf("./log/compute_dist/%s/%s", opt$proj_name, opt$task_id), recursive=TRUE)
saveRDS(opt, sprintf("./log/compute_dist/%s/%s/opt_handles.rds", opt$proj_name, opt$task_id))

## Read in data
source(file.path("scripts", "utils", "utils.R"))

gene_exp_train <- try(read.table(opt$gene_train, sep="\t",header=TRUE))
if (class(gene_exp_train) == "try-error"){
  message("Training gene expression file not found: Please input a valid training gene expression file location! ")
  gene_exp_train
}
if(opt$verbose){
  message("Data: gene_exp_train")
}

colnames(gene_exp_train) <- unlist(lapply(colnames(gene_exp_train), function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x = str_replace_all(x, "\\+", ".")
  x
}))

gene_exp_test <- try(read.table(opt$gene_test, sep="\t",header=TRUE))
if (class(gene_exp_test) == "try-error"){
  message("Testing gene expression file not found: Please input a valid testing gene expression file location! ")
  gene_exp_test
}
if(opt$verbose){
  message("Data: gene_exp_test")
}
colnames(gene_exp_test) <- unlist(lapply(colnames(gene_exp_test), function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x = str_replace_all(x, "\\+", ".")
  x
}))

test_only_samples <- setdiff(colnames(gene_exp_test), colnames(gene_exp_train))
if(opt$verbose){
  message(sprintf("Number of samples not in reference: %s", length(test_only_samples)))
}
if (length(test_only_samples) > 0){
  train_genes <- row.names(gene_exp_train)
  test_genes <- row.names(gene_exp_test)
  common_genes <- intersect(train_genes, test_genes)
  if (length(common_genes) == 0){
    message("No common genes found between training and testing datasets. Please check your input files.")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  gene_exp_train <- gene_exp_train[common_genes, , drop=FALSE]
  gene_exp_test <- gene_exp_test[common_genes, , drop=FALSE]

  gene_exp <- cbind(gene_exp_train, gene_exp_test[,test_only_samples])
  gene_exp_mat <- as.matrix(gene_exp)
  normalized_data <- normalize.quantiles(gene_exp_mat)
  gene_exp_qn <- as.data.frame(normalized_data)
  row.names(gene_exp_qn) <- row.names(gene_exp)
  colnames(gene_exp_qn) <- colnames(gene_exp)
} else {
  gene_exp_qn <- gene_exp_train
}

proj_name <- opt$proj_name
out_dir <- opt$outdir
dist_method <- opt$dist

gene_exp_mat_qn <- as.matrix(gene_exp_qn)
gene_exp_mat_qn_st <- scale(gene_exp_mat_qn)

## Using all gene expression for measuring distance
dist_save <- file.path(out_dir, proj_name, "dist_nearest_ct")
dir.create(dist_save, recursive = TRUE)

#### Compute Euclidean Distance
if(dist_method == "euc"){
  if(!file.exists(file.path(dist_save, "dist_ct_df.rds"))){
    dist_mat = compute_euc_distance(gene_exp_mat_qn_st, gene_exp_mat_qn_st)
    saveRDS(dist_mat, file.path(dist_save, "dist_ct_df.rds"))
  }
  if(opt$verbose){
    message("Euclidean distance matrix computed.")
  }
  dist_mat <- readRDS(file.path(dist_save, "dist_ct_df.rds"))
}

#### Compute Correlation
if(dist_method != "euc"){
  if(dist_method != "scc" & dist_method != "pcc"){
    message("Please specify a distance method of one of the following: 'euc', 'pcc' or 'scc'!")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  if(!file.exists(file.path(dist_save, sprintf("corr_ct_df.rds", dist_method)))){
    corr_mat = compute_corr(gene_exp_mat_qn_st, gene_exp_mat_qn_st, metric=dist_method)
    saveRDS(corr_mat, file.path(dist_save, sprintf("corr_ct_df.rds", dist_method)))
  }
  if(opt$verbose){
    message("Correlation distance matrix computed")
  }
  corr_mat <- readRDS(file.path(dist_save, sprintf("corr_ct_df.rds", dist_method)))
}

#### Get nearest celltypes using one of the metrics
all_ct <- unique(c(colnames(gene_exp_train), colnames(gene_exp_test)))

if(!file.exists(file.path(dist_save, sprintf("all_gene_dict.rds", dist_method)))){
  all_gene_dict <- list()
  for(ct in all_ct){
    if(dist_method == "euc"){
      sorted_ct = names(sort(dist_mat[,ct]))
    } else {
      sorted_ct = names(sort(corr_mat[,ct], decreasing = TRUE, na.last = TRUE))
    }
    all_gene_dict[[ct]] <- sorted_ct
  }

  saveRDS(all_gene_dict, file.path(dist_save, sprintf("all_gene_dict.rds", dist_method)))
}
if(opt$verbose){
  message("Sample distances calculated.")
}