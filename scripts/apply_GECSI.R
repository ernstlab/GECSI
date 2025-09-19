#!/usr/bin/env Rscript
### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
options(scipen=999) ## avoid scientific notation
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# getwd()
gc()
memory.limit(size=Inf)

#.libPaths(c("/Library/Frameworks/R.framework/Versions/3.6/Resources/library", .libPaths()))
#print(.libPaths())
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, ggplot2, GenomicRanges, caret, foreach, doParallel, data.table, rtracklayer, glmnet, pryr, fastDummies) # parallel

parser <- OptionParser()
option_list <- list(
  make_option(c("-v", "--verbose"), type="logical", default=FALSE,
              help="Print extra output"),
  make_option(c("--ignore_missing_model"), type="logical", default=FALSE, help="If a model for a sample and state is not found, ignore it and predict all zeros. Otherwise skip the sample and go to next sample."),
  make_option(c("--overwrite"), type="logical", default=FALSE, help="If true, then overwrite all existing files in the project folder."),
  make_option(c("--train_chr"), default="all",
              help = "Which chromosome to use (\"chr1\"-\"chrX\", or \"all\" if use all chromosomes) [default \"%default\"]",
              metavar="train chromosome"),
  make_option(c("--test_chr"), default="chr1",
              help = "Which chromosome to use (\"chr1\"-\"chrX\")",
              metavar="test chromosome"),
  make_option(c("--num_train_ct"), type="integer", default=40, help="Number of training cell type models to use for prediction [default %default]. If num_train_ct='all', training_samples parameter need to be specified.", metavar="num train celltypes"),
  make_option(c("--holdout"), default="None",
              help = "Held out cell type (or an rds file specifying a list of held out cell types). If held out cell type in the list of training celltypes, remove it from training.",
              metavar="holdout"),
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c("-c", "--chip"), default="./data/chip_dir/",
              help = "Set directory for ChIP-seq data [default \"%default\"]",
              metavar="chip-seq data directory"),
  make_option(c("--proj_name"), default="new_proj",
              help = "Define a project name. Distance metrics will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name"),
  make_option(c("--parallel"), type="integer", default=3, help = "Number of cores for parallel processing reference data.", metavar="parallel"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--seed"), default=42,
              help = "Set seed for this experiment (such as training-validation split)",
              metavar="seed"),
  make_option(c("--sample_size"), default=500000, help = "Sample size for data reading into the model", metavar="sample size"),
  make_option(c("--num_states"), default=18, help = "Number of total chromatin states the model is trained to predict", metavar="number of chromatin states"),
  make_option(c("--states_list"), default="", help = "A file providing a list of available chromatin state names in order, separated in lines. If not provided, will default to the Roadmap 18-state chromatin state ordering.", metavar="States list"),
  make_option(c("--task_id"), type="integer", default=1, help = "Task ID for this experiment", metavar = "task id"),
  make_option(c("--dist"), default="pcc", help = "Distance method", metavar="distance method"),
  make_option(c("-k", "--num_nearest_celltypes"), type="integer", default=5, help="Number of k nearest training cell types to look for [default %default]", metavar="num nearest celltypes"),
  make_option(c("--training_samples"), default="./training_samples", help="Path to training samples", metavar="path to training samples"),
  make_option(c("--lambda"), default=0.001, help="Lambda for Lasso regularization", metavar="lasso regularization parameter"),
  make_option(c("--probability"), type="logical", default=FALSE, help="Generate probability estimation files", metavar="generate probability estimation")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(opt$verbose){
  print(opt)
}

dir.create(sprintf("./log/apply/%s/%s", opt$proj_name, opt$task_id), recursive=TRUE)
saveRDS(opt, sprintf("./log/apply/%s/%s/opt_handles.rds", opt$proj_name, opt$task_id))

### TAKE INPUT FROM ARGUMENT HANDLES
wd <- opt$wd
setwd(wd)
source(file.path("scripts", "utils", "utils.R"))

train_chr <- opt$train_chr
test_chr <- opt$test_chr
holdout_f <- opt$holdout
k <- opt$num_nearest_celltypes
dist_method <- opt$dist
lambda <- opt$lambda
nstates <- opt$num_states
states_list_f <- opt$states_list
generate_probability <- opt$probability

if (holdout_f != "None"){
  if(!endsWith(holdout_f, ".txt")){
    holdout <- holdout_f
    holdout = str_replace_all(holdout, " ", ".")
    holdout = str_replace_all(holdout, "-", ".")
    holdout = str_replace_all(holdout, "\\+", ".")
  } else {
    holdout <- try(readLines(holdout_f))
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

if (states_list_f != ""){
  states_list <- try(readLines(states_list_f))
  if (class(states_list) == "try-error"){
    message("States list file not found or not in a compatible format!")
    states_list
  }
  states_list <- as.character(states_list)
  states_list <- unlist(lapply(states_list, function(x){
    x = str_replace_all(x, " ", ".")
    x = str_replace_all(x, "-", ".")
    x = str_replace_all(x, "\\+", ".")
    x
  }))
} else {
  states_list <- c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk","7_EnhG1","8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF/Rpts","13_Het","14_TssBiv","15_EnhBiv","16_ReprPC","17_ReprPCWk","18_Quies")
}

output_dir <- opt$outdir
proj_name <- opt$proj_name
num_train_ct <- opt$num_train_ct
sample_size <- opt$sample_size
task_id <- opt$task_id

proj_name_train <- sprintf("%s/Train_%s", proj_name, task_id)
proj_name_test <- sprintf("%s/Apply_%s", proj_name, task_id)

dir.create(file.path(output_dir, proj_name_test), recursive = TRUE)

seed <- opt$seed
pc <- opt$pc

### Read training cell types for this PC and this holdout sample
train_nearest_ct_save <- file.path(output_dir, proj_name, "find_nearest_ct", sprintf("nearest_ct_%s", task_id))
if (num_train_ct == "all"){
  training_ct <- readLines(opt$training_samples)
} else {
  training_ct <- readLines(file.path(train_nearest_ct_save, sprintf("samples_to_train_%s_%s_%s.txt", num_train_ct, dist_method, holdout)))
}

all_training_samples <- readLines(opt$training_samples)

dist_save <- file.path(output_dir, proj_name, "dist_nearest_ct")
dist_dict <- readRDS(file.path(dist_save, sprintf("all_gene_dict.rds", dist_method)))
dist_ss <- dist_dict
dist_ct <- dist_ss[[holdout]]
dist_ct_use <- dist_ct[which(dist_ct %in% all_training_samples)]
dist_ct_rm <- setdiff(dist_ct_use, holdout)

k_training_ct <- dist_ct_rm

### READ IN CHROMATIN STATE FILES FOR TRAINING CELL TYPES
## Generate GR first
## Generate a ChIP-seq chunk from a list of samples first
all_chip_files <- list.files(opt$chip, full.names = TRUE)
ct = training_ct[1]
pattern <- paste0(ct, "([_\\-\\. ]|$)")
read_chip = grep(pattern, all_chip_files, value = TRUE)
if(length(read_chip) == 0){
  stop(paste0("Correponding ChIP-seq file not found for cell type ", ct, "!"))
} else if(length(read_chip) > 1){
  stop(paste0("Multiple ChIP-seq files found for cell type ", ct, ": ", paste(read_chip, collapse=", ")))
}
if(endsWith(read_chip, ".rds")){
  chip_ct = readRDS(read_chip)
} else if(endsWith(read_chip, ".bed") || endsWith(read_chip, ".bed.gz")){
  chip_ct = import(read_chip, format="bed")
  names(mcols(chip_ct)) <- ct
} else {
  message("ChIP-seq file not in a compatible format!")
  chip_ct
}

if (test_chr != "all"){
  chip_ct = chip_ct[seqnames(chip_ct) == test_chr,]
}

remove_slash <- function(x){
  if(grepl("/", x)){
    x <- strsplit(x, split="/")[[1]][1]
  }
  return(x)
}

k_list <- paste0("k", seq(1,k))

chip_test_path <- file.path(output_dir, proj_name_test, "knn-tracks-test", test_chr)
chip_test_save <- sprintf("ATLAS_CHIP_test_%s_%s.rds", holdout, k)

if (!file.exists(file.path(chip_test_path, chip_test_save)) | opt$overwrite){
    chip_ct_chunk = chip_ct
    chip_ct_chunk_gr = granges(chip_ct_chunk)

    no_cores <- opt$parallel # number of parallel cores
    start.t.total<-Sys.time()
    if(opt$verbose){
      print(sprintf("Starting parallel processing with %s cores", no_cores))
      print(sprintf("Start time: %s", start.t.total))
      print("Generating chromatin state features for heldout cell type...")
    }
    startt<-Sys.time()

    cl<-makeCluster(no_cores,outfile = file.path(output_dir, proj_name_test, sprintf("myCalculationDebug_y_train_%s_%s.txt", holdout, pc)))
    registerDoParallel(cl)

    generate_chunk <- function(ct){
        all_chip_files <- list.files(opt$chip, full.names = TRUE)
        pattern <- paste0(ct, "([_\\-\\. ]|$)")
        read_chip = grep(pattern, all_chip_files, value = TRUE)
        if(length(read_chip) == 0){
          stop(paste0("Correponding ChIP-seq file not found for cell type ", ct, "!"))
        } else if(length(read_chip) > 1){
          stop(paste0("Multiple ChIP-seq files found for cell type ", ct, ": ", paste(read_chip, collapse=", ")))
        }
        if(endsWith(read_chip, ".rds")){
          chip_ct = readRDS(read_chip)
        } else if(endsWith(read_chip, ".bed") || endsWith(read_chip, ".bed.gz")){
          chip_ct = import(read_chip, format="bed")
        } else {
          message("ChIP-seq file not in a compatible format!")
          chip_ct
        }
        if(opt$verbose){
          print(ct)
          print(read_chip)
        }
        if(test_chr != "all"){
          chip_ct = chip_ct[seqnames(chip_ct) == test_chr,]
        }
        names(mcols(chip_ct)) <- ct
        chip_ct_chunk = chip_ct
        return(as.data.frame(mcols(chip_ct_chunk[,ct])))
    }

    ### TODO
    if(k < 15){
      chip_ct_chunk_df <- foreach(train_ct=1:15
                          ,.combine=cbind # data from different threads will be in one table
                          ,.packages=c("GenomicRanges", "rtracklayer")# All packages that your funtion is using must be called here
                          ,.inorder=T) %dopar% #don`t forget this directive
                          {
                          generate_chunk(k_training_ct[train_ct])
                          }
      stopCluster(cl)
    } else {
      chip_ct_chunk_df <- foreach(train_ct=1:k
                          ,.combine=cbind # data from different threads will be in one table
                          ,.packages=c("GenomicRanges", "rtracklayer")# All packages that your funtion is using must be called here
                          ,.inorder=T) %dopar% #don`t forget this directive
                          {
                          generate_chunk(k_training_ct[train_ct])
                          }
      stopCluster(cl)
    }


    # chip_ct_chunk_df <- do.call(cbind, chip_ct_chunk_l)
    mcols(chip_ct_chunk_gr) <- chip_ct_chunk_df

    ATLAS_CHIP <- chip_ct_chunk_gr
    X_test <- ATLAS_CHIP[, k_training_ct[1:k], drop = FALSE]
    dir.create(chip_test_path, recursive=TRUE)

    if (k < 15){
      X_test_save <- ATLAS_CHIP[, k_training_ct[1:15], drop = FALSE]
    } else {
      X_test_save <- X_test
    }
    saveRDS(X_test_save, file.path(chip_test_path, chip_test_save))

}
if(opt$verbose){
  print(sprintf("Features for heldout %s saved!", holdout))
}
rm(chip_ct)
rm(ATLAS_CHIP)
gc()

X_test <- readRDS(file.path(chip_test_path, chip_test_save))
# ATLAS_CHIP <- X_test

X_test_df <- as.data.frame(mcols(X_test))
if(k < 15){
  X_test_df <- X_test_df[,1:k]
}
names(X_test_df) <- k_list

# test_1hot_df <- as.data.frame(model.matrix(~ . - 1, data = X_test_df))
test_1hot_df <- dummy_cols(X_test_df, remove_first_dummy = FALSE, remove_selected_columns = TRUE)

colnames(test_1hot_df) <- sub("^(k\\d+)_", "\\1", colnames(test_1hot_df))
colnames(test_1hot_df) <- unlist(lapply(colnames(test_1hot_df), remove_slash))

if(opt$verbose){
  print("One-hot encoding generated!")
}

all_states <- sort(unique(X_test_df[,1]))
all_states <- unlist(lapply(all_states, remove_slash))
nstate = length(all_states)
combinations <- expand.grid(k_list = k_list, all_states = all_states)
all_predictors <- paste0(combinations$k_list, combinations$all_states)
formula_knn <- paste(all_predictors, collapse = ' + ')

predictors <- intersect(colnames(test_1hot_df), all_predictors)
pred_dir <- file.path(output_dir, proj_name_test, "predictions", "lr-model", test_chr)
dir.create(pred_dir, recursive=TRUE)

all_states = states_list

all_ct_prob_df = data.frame(matrix(0,nrow=nrow(test_1hot_df), ncol=nstates))
names(all_ct_prob_df) = all_states

states_prob_dir = file.path(pred_dir, sprintf("states_prob_df_combined_train%s_k%s_s%s_%s_lambda%s.rds", num_train_ct, k, sample_size, holdout, lambda))

if((!file.exists(states_prob_dir)) | opt$overwrite){
    get_state_probability <- function(train_ct) {
      model_save = file.path(output_dir, proj_name_train, "models", "lr-model-multinom", train_chr, sprintf("lr_model_k%s_global_val_%s_s%s.rds", k, train_ct, sample_size))
      model <- try(readRDS(model_save))
      if(class(model) == "try-error"){
        message(sprintf("Model file for training celltype %s for PC %s not found!! ", train_ct, pc))
        if (opt$ignore_missing_model){
            return(0)            
        } else {
            stop("Error returned. If you want to ignore missing model, set option handle --ignore_missing_model")
        }            
      }

      coefn <- rownames(coef(model)$X1)[-1]
      add_coef <- setdiff(coefn, colnames(test_1hot_df))
      test_1hot_df[add_coef] <- 0
      remove_coef <- setdiff(colnames(test_1hot_df), coefn)
      newx <- test_1hot_df[,coefn]

      pp_test <- predict(model, newx=as.matrix(newx), s=lambda, type = "response")
      new_colnames <- factor(colnames(pp_test), levels=paste0("X",seq(1,length(all_states))))
      levels(new_colnames) <- all_states
      colnames(pp_test) <- as.character(new_colnames)
      if(opt$verbose){
        print(sprintf("Prediction from %s generated!", train_ct))
      }
      return(pp_test)
    }

    for(train_ct in training_ct){
      prob_df <- get_state_probability(train_ct)
      if(!identical(prob_df, 0)){
        # print(prob_df)
        prob_df <- as.data.frame(prob_df[,,1])
        add_names <- setdiff(all_states, colnames(prob_df))
        prob_df[add_names] <- 0
        prob_df <- prob_df[,all_states]
        all_ct_prob_df = all_ct_prob_df + prob_df
      }
    }
    saveRDS(all_ct_prob_df, states_prob_dir)
}

rm(test_1hot_df)
# rm(X_test)
rm(X_test_df)
gc()
states_prob_df <- readRDS(states_prob_dir)

if(opt$probability){
  prob_dir <- file.path(pred_dir, "probability_estimate_files")
  dir.create(prob_dir, recursive = TRUE)
  if(test_chr == "all"){
    for(chr in paste0("chr", c(seq(1,22), "X"))){
      print(chr)
      output_file = file.path(prob_dir, sprintf("%s_%s_%s_probability_estimate.txt.gz", holdout, length(all_states), chr))
      states_prob_df = states_prob_df / length(training_ct)

      # Write to output file
      conn <- gzfile(output_file, "w")

      # Write first line: region ID and chromosome
      writeLines(paste(holdout, chr, sep = "\t"), conn)

      # Write second line: header (E1, E2, ..., E18)
      writeLines(paste(colnames(states_prob_df), collapse = "\t"), conn)

      # Write posterior probabilities
      apply(states_prob_df, 1, function(row) {
          writeLines(paste(row, collapse = "\t"), conn)
      })

      # Close connection
      close(conn)
    }
  }
  else {
    output_file = file.path(prob_dir, sprintf("%s_%s_%s_probability_estimate.txt.gz", holdout, length(all_states), test_chr))
    states_prob_df = states_prob_df / length(training_ct)

    # Write to output file
    conn <- gzfile(output_file, "w")

    # Write first line: region ID and chromosome
    writeLines(paste(holdout, test_chr, sep = "\t"), conn)

    # Write second line: header (E1, E2, ..., E18)
    writeLines(paste(colnames(states_prob_df), collapse = "\t"), conn)

    # Write posterior probabilities
    apply(states_prob_df, 1, function(row) {
        writeLines(paste(row, collapse = "\t"), conn)
    })

    # Close connection
    close(conn)
  }
  if (opt$verbose){
  print("Probability estimation files generated!")
  }
}

pred_states_dir <- file.path(pred_dir, sprintf("pred_states_train%s_k%s_s%s_%s_lambda%s.rds", num_train_ct, k, sample_size, holdout, lambda))
if((!file.exists(pred_states_dir)) | opt$overwrite){
    highest_prob_state_in_row <- function(row) {
      most_frequent_state <- names(row[which.max(row)])
      return(most_frequent_state)
    }

    most_frequent_states <- lapply(1:nrow(states_prob_df), function(row){
                                if(row %% 1000 == 0) { if (opt$verbose) {print(row)}}
                                highest_prob_state_in_row(states_prob_df[row,])
                            })
    # most_frequent_states <- apply(states_prob_df, 1, most_frequent_state_in_row)
    saveRDS(most_frequent_states, pred_states_dir)
}
if(opt$verbose){
  print("Most frequent states generated!")
}

pred_bed_dir <- file.path(pred_dir, "bed_files", sprintf("%s_pred_states_train%s_k%s_s%s_lambda%s.bed.gz", holdout, num_train_ct, k, sample_size, lambda))
if((!file.exists(pred_bed_dir)) | opt$overwrite){
  most_frequent_states <- readRDS(pred_states_dir)
  most_frequent_states <- unlist(most_frequent_states)
  # gr <- readRDS("../data/granges_200bp_bins.rds")
  gr <- granges(X_test)
  if(test_chr != "all"){
    gr <- gr[seqnames(gr) == test_chr,]
  }
  mcols(gr) = most_frequent_states
  
  names(mcols(gr)) <- "State"
  unique_states <- unlist(reduce(split(gr, ~State)))
  unique_states$State = names(unique_states)
  unique_states <- sortSeqlevels(unique_states)
  unique_states <- sort(unique_states)
  names(unique_states) <- NULL

  most_frequent_states = unique_states$State
  gr <- unique_states

  bed_df <- data.frame(matrix(nrow=length(most_frequent_states), ncol=8))
  names(bed_df) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd")

  bed_df$chrom = as.factor(seqnames(gr))
  bed_df$chromStart = start(gr)-1
  bed_df$chromEnd = end(gr)
  bed_df$name = most_frequent_states
  bed_df$score = 0
  bed_df$strand = "."
  bed_df$thickStart = bed_df$chromStart
  bed_df$thickEnd = bed_df$chromEnd
  color_scheme = read.table("./data/color_scheme.txt", sep= "\t", header=TRUE)
  color_scheme$name = paste0(color_scheme$`STATE.NO.`, "_", color_scheme$`MNEMONIC`)
  color_scheme$itemRgb = color_scheme$COLOR.CODE
  color_scheme <- color_scheme[,c("name","itemRgb")]
  bed_df <- merge(bed_df, color_scheme, by="name")
  bed_df <- bed_df[,c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")]
  bed_df <- bed_df[order(bed_df$chromStart), ]

  bed_dir <- file.path(pred_dir, "bed_files")
  dir.create(bed_dir, recursive = TRUE)
  write.table(bed_df, gzfile(pred_bed_dir), sep = "\t", quote=FALSE, col.names = FALSE, row.names=FALSE)
}

