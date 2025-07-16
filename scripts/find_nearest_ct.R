#!/usr/bin/env Rscript
### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# getwd()
gc()
memory.limit(size=Inf)

#.libPaths(c("/Library/Frameworks/R.framework/Versions/3.6/Resources/library", .libPaths()))
#print(.libPaths())
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, data.table) # parallel

parser <- OptionParser()
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output [default]"),
  make_option(c("--training_ct"), default="./path_to_training_ct/training_ct.txt",
              help = "A file listing celltypes to use for training [default \"%default\"]",
              metavar="training celltypes"),
  make_option(c("--holdout"), default="None",
              help = "Held out cell types (an rds file specifying a list of held out cell types). If held out cell types are in the list of training celltypes, remove it from training.",
              metavar="holdout"),
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c("--proj_name"), default="new_proj",
              help = "The project name folder where the distance metrics are stored, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--num_ct_to_train"), default=40,
              help = "Number of nearest ct to train",
              metavar="num_ct_to_train"),
  make_option(c("--dist"), default="pcc",
            help = "Method to use to calculate nearest training cell-types. 'pcc' for Pearson Correlation, 'scc' for Spearman Correlation, and 'euc' for Euclidean distance)",
            metavar="dist"),
  make_option(c("--task_id"), default="1", help = "Task identifier (user-specified) for this experiment", metavar = "task id")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(opt$verbose){
  print(opt)
}

dir.create(sprintf("./log/find_nearest_ct/%s/%s", opt$proj_name, opt$task_id), recursive=TRUE)
saveRDS(opt, sprintf("./log/find_nearest_ct/%s/%s/opt_handles.rds", opt$proj_name, opt$task_id))


### Get parameters
wd <- opt$wd
setwd(wd)
source(file.path("scripts", "utils", "utils.R"))


training_ct <- try(readLines(opt$training_ct))
if (class(training_ct) == "try-error"){
  message("Training celltypes file not found!!")
  training_ct
}
training_ct <- as.character(training_ct)
training_ct <- unlist(lapply(training_ct, function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x = str_replace_all(x, "\\+", ".")
  x
}))

num_train_ct <- length(training_ct)

proj_name <- opt$proj_name
across_dir <- opt$outdir
dist_method <- opt$dist
holdout_f <- opt$holdout
dist_method <- opt$dist
task_id <- opt$task_id


num_ct_to_train <- opt$num_ct_to_train

if (holdout_f != "None"){
  if(!endsWith(holdout_f, ".txt")){
    holdout <- holdout_f
    holdout = str_replace_all(holdout, " ", ".")
    holdout = str_replace_all(holdout, "-", ".")
    holdout = str_replace_all(holdout, "\\+", ".")
    holdout_id = holdout
  } else {
    holdout <- try(readLines(holdout_f))
    # holdout_id <- tools::file_path_sans_ext(basename(holdout_f))
    if (class(holdout) == "try-error"){
      message("Held-out celltypes file not found!!")
      holdout
    }
    holdout <- as.character(holdout)
    holdout <- unlist(lapply(holdout, function(x){
      x = str_replace_all(x, " ", ".")
      x = str_replace_all(x, "-", ".")
      x = str_replace_all(x, "\\+", ".")
      x
    }))
  }
}

## Read PCA-based and overall gene expression-based nearest cell types

dist_save <- file.path(across_dir, proj_name, "dist_nearest_ct")
dist_dict <- try(readRDS(file.path(dist_save, sprintf("all_gene_dict.rds"))))
if (class(dist_dict) == "try-error"){
  message("PCA-based nearest sample file not found!!")
  dist_dict
}
dists <- names(dist_dict)

train_nearest_ct_save <- file.path(across_dir, proj_name, "find_nearest_ct", sprintf("nearest_ct_%s",task_id))
dir.create(train_nearest_ct_save, recursive=TRUE)

num_test_ct <- length(holdout)

for (test_ct in holdout){
      if(!file.exists(file.path(train_nearest_ct_save, sprintf("samples_to_train_%s_%s_%s.txt", num_ct_to_train, dist_method, test_ct)))){
        dist_ct <- dist_dict[[test_ct]]
        dist_ct <- dist_ct[which(dist_ct %in% c(training_ct))]
        dist_ct_train <- dist_ct[which(!dist_ct %in% c(holdout))][1:num_ct_to_train]
        writeLines(dist_ct_train, file.path(train_nearest_ct_save, sprintf("samples_to_train_%s_%s_%s.txt", num_ct_to_train, dist_method, test_ct)))
        # writeLines(remove_version(dist_ct_train), file.path(train_nearest_ct_save, sprintf("samples_to_train_%s_%s_holdout_%s_fold%s.txt", num_ct_to_train, dist_method, test_ct, fold)))
    }
}

if (opt$verbose){
  message("Nearest samples are found and saved!")
}