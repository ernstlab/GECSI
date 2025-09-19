#!/usr/bin/env Rscript
### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
options(scipen = 999)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
# getwd()
gc()
memory.limit(size=Inf)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))


### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, ggplot2, GenomicRanges, caret, foreach, doParallel, data.table, rtracklayer, glmnet) # parallel

parser <- OptionParser()
option_list <- list(
  make_option(c("-v", "--verbose"), type="logical", default=FALSE,
              help="Print extra output"),
  make_option(c("--overwrite"), type="logical", default=FALSE, help="If true, then overwrite all existing files in the project folder."),
  make_option(c("--chr"), default="all",
              help = "Which chromosome to use (\"chr1\"-\"chrX\", or \"all\" if use all chromosomes) [default \"%default\"]",
              metavar="chromosome"),
  make_option(c("--training_ct"), default="./path_to_training_ct/training_ct.txt",
              help = "Training cell type (or an rds file listing celltypes to use for training [default \"%default\"])",
              metavar="training celltype"),
  make_option(c("--all_training_samples"), default="None",
              help = "All training samples that can be used for training this cell type (an rds file specifying a list of training cell types). ",
              metavar="all train"),
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c("-c", "--chip"), default="./data/chip_dir/",
              help = "Set directory for ChIP-seq data [default \"%default\"]",
              metavar="chip-seq data directory"),
  make_option(c("--proj_name"), default="new_proj",
              help = "Define a project name. Distance metrics will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--seed"), default=42,
              help = "Set seed for this experiment (such as training-validation split)",
              metavar="seed"),
  make_option(c("--sample_size"), default=100000, help = "Sample size for data reading into the model", metavar="sample size"),
  make_option(c("--num_states"), default=18, help = "Number of total chromatin states the model is trained to predict", metavar="number of chromatin states"),
  make_option(c("--quies_state"), default="18_Quies", help = "The name of the quiescent chromatin state. For example in Roadmap 18 state model, this state is called '18_Quies'", metavar="Name of quiescent state"),
  make_option(c("--states_list"), default="", help = "A file providing a list of available chromatin state names in order, separated in lines. If not provided, will default to the Roadmap 18-state chromatin state ordering.", metavar="States list"),
  make_option(c("--parallel"), type="integer", default=3, help = "Number of cores for parallel processing reference data.", metavar="parallel"),
  make_option(c("--task_id"), default="1", help = "Task identifier (user-specified) for this experiment", metavar = "task id"),
  make_option(c("--dist"), default="pcc", help = "Distance method", metavar="distance method"),
  make_option(c("-k", "--num_nearest_celltypes"), type="integer", default=10, help="Number of k nearest training cell types to look for [default %default]", metavar="num nearest celltypes"),
  make_option(c("--save_train"), type="logical", default=FALSE,
              help = "Use this flag if you want to save training predictions. Otherwise only the models will be saved.")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(opt$verbose){
  print(opt)
}

dir.create(sprintf("./log/train/%s/%s", opt$proj_name, opt$task_id), recursive=TRUE)
saveRDS(opt, sprintf("./log/train/%s/%s/opt_handles.rds", opt$proj_name, opt$task_id))

### TAKE INPUT FROM ARGUMENT HANDLES
wd <- opt$wd
setwd(wd)
source(file.path("scripts", "utils", "utils.R"))

chr <- opt$chr
training_ct_f <- opt$training_ct
all_training_samples_f <- opt$all_training_samples
k <- opt$num_nearest_celltypes
dist_method <- opt$dist
task_id <- opt$task_id
nstates <- opt$num_states
quies_state <- opt$quies_state
states_list_f <- opt$states_list

if (training_ct_f != "None"){
  if(!endsWith(training_ct_f, ".txt")){
    training_ct <- training_ct_f
    training_ct = str_replace_all(training_ct, " ", ".")
    training_ct = str_replace_all(training_ct, "-", ".")
    training_ct = str_replace_all(training_ct, "\\+", ".")
  } else {
    training_ct <- try(readLines(training_ct_f))
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
  }
}

if (all_training_samples_f != "None"){
    all_training_samples <- try(readLines(all_training_samples_f))
    if (class(all_training_samples) == "try-error"){
        message("A list of reference samples file not found or not in a compatible format!")
        all_training_samples
    }
    all_training_samples <- as.character(all_training_samples)
    all_training_samples <- unlist(lapply(all_training_samples, function(x){
        x = str_replace_all(x, " ", ".")
        x = str_replace_all(x, "-", ".")
        x = str_replace_all(x, "\\+", ".")
        x
    }))
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

num_train_ct <- length(training_ct)

output_dir <- opt$outdir
proj_name <- opt$proj_name
proj_name_train <- sprintf("%s/Train_%s", proj_name, task_id)
dir.create(file.path(output_dir, proj_name_train), recursive = TRUE)

seed <- opt$seed

### SET WORKING DIRECTORY AND FORMAT INPUT DATA
chip_path = opt$chip
sample_size = as.numeric(opt$sample_size)

### READ IN CHIP SEQ DATA FROM FOLDER

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

if (chr != "all"){
  chip_ct = chip_ct[seqnames(chip_ct) == chr,]
}


## Now generate chip seq data for a set of sampled positions
if (!file.exists(file.path(output_dir, proj_name, "sampled_pos", sprintf("%s_%s_sampled_pos_%s.rds", sample_size, chr, training_ct)))){
  set.seed(seed)
  all_pos <- unname(unlist(mcols(chip_ct)))
  state_freq <- table(all_pos)
  sorted_table <- sort(state_freq, decreasing=FALSE)
  sampled_pos = c()
  num_states = length(names(sorted_table))
  for(state_name in names(sorted_table)){
    if(state_name == names(sorted_table)[num_states]){
      num_add_pos = sample_size - length(sampled_pos)
      add_pos = sample(which(all_pos == state_name), num_add_pos)
      sampled_pos <- c(sampled_pos, add_pos)
    } else {
      num_add_pos = state_freq[state_name]/sum(state_freq)*sample_size
      add_pos = sample(which(all_pos == state_name), num_add_pos)
      sampled_pos <- c(sampled_pos, add_pos)
    }
  }
  dir.create(file.path(output_dir, proj_name, "sampled_pos"), recursive=TRUE)
  saveRDS(sampled_pos, file.path(output_dir, proj_name, "sampled_pos", sprintf("%s_%s_sampled_pos_%s.rds", sample_size, chr, training_ct)))
}


sampled_pos <- readRDS(file.path(output_dir, proj_name, "sampled_pos", sprintf("%s_%s_sampled_pos_%s.rds", sample_size, chr, training_ct)))

dist_save <- file.path(output_dir, proj_name, "dist_nearest_ct")
dist_dict <- readRDS(file.path(dist_save, sprintf("all_gene_dict.rds", dist_method)))
dist_ss <- dist_dict

remove_slash <- function(x){
  if(grepl("/", x)){
    x <- strsplit(x, split="/")[[1]][1]
  }
  return(x)
}

k_list <- paste0("k", seq(1,k))

for (ct in training_ct){
  dist_ct <- dist_ss[[ct]]
  dist_ct_use <- intersect(dist_ct, all_training_samples)
  dist_ct_rm <- setdiff(dist_ct_use, ct)
  dist_ct_train <- dist_ct_rm[1:k]
  dist_ct_train_w_tg <- c(ct, dist_ct_train)

  chip_train_path <- file.path(output_dir, proj_name_train, "knn-tracks-train", chr)
  chip_train_save <- sprintf("ATLAS_CHIP_train_%s_%s_%s.rds", sample_size, ct, k)
  if ((!file.exists(file.path(chip_train_path, chip_train_save))) | opt$overwrite){
    chip_ct_chunk = chip_ct[sampled_pos,]
    chip_ct_chunk_gr = granges(chip_ct_chunk)

    no_cores <- opt$parallel
    if(opt$verbose){
      start.t.total<-Sys.time()
      print(start.t.total)
      startt<-Sys.time()
      print(startt)
    }

    cl<-makeCluster(no_cores,outfile = file.path(output_dir, proj_name_train, sprintf("myCalculationDebug_y_train_%s.txt", ct)))
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
      if (chr != "all"){
        chip_ct = chip_ct[seqnames(chip_ct) == chr,]
      }
      names(mcols(chip_ct)) <- ct
      chip_ct_chunk = chip_ct[sampled_pos,]
      return(as.data.frame(mcols(chip_ct_chunk[,ct])))
    }

    ### TODO
    chip_ct_chunk_df <- foreach(train_ct=1:(k+1)
                      ,.combine=cbind # data from different threads will be in one table
                      ,.packages=c("GenomicRanges", "rtracklayer")# All packages that your funtion is using must be called here
                      ,.inorder=T) %dopar% #don`t forget this directive
                      {
                        generate_chunk(dist_ct_train_w_tg[train_ct])
                      }
    
    stopCluster(cl)

    mcols(chip_ct_chunk_gr) <- chip_ct_chunk_df

    ATLAS_CHIP <- chip_ct_chunk_gr
    Y_train <- ATLAS_CHIP[, dist_ct_train_w_tg, drop = FALSE]

    dir.create(chip_train_path, recursive=TRUE)
    saveRDS(Y_train, file.path(chip_train_path, chip_train_save))

  }

  if(opt$verbose){
    print("Features saved!")
  }

  Y_train <- readRDS(file.path(chip_train_path, chip_train_save))
  ATLAS_CHIP <- Y_train

  Y_train_df <- as.data.frame(mcols(Y_train))
  names(Y_train_df) <- c(ct, k_list)
  Y_train_df_predictor <- Y_train_df[,-1]

  train_1hot_df <- as.data.frame(model.matrix(~ . - 1, data = Y_train_df_predictor))
  colnames(train_1hot_df) <- unlist(lapply(colnames(train_1hot_df), remove_slash))
  train_1hot_df <- train_1hot_df[,names(which(colSums(train_1hot_df) > 10))]
  if(opt$verbose){
    print("One-hot encoding generated!")
  }

  target <- ct
  predictors <- colnames(train_1hot_df)

  train_1hot_df[,target] <- Y_train_df[,target]
  target_factor <- as.factor(train_1hot_df[,target])
  if (quies_state %in% levels(target_factor)) {
    train_1hot_df[,target] <- relevel(target_factor, ref = quies_state)
  } else {
    # Find the most frequent state
    most_freq <- names(sort(table(target_factor), decreasing = TRUE))[1]
    message(paste0("quies_state '", quies_state, "' not found. Using most frequent state '", most_freq, "' as reference."))
    train_1hot_df[,target] <- relevel(target_factor, ref = most_freq)
  }

  formula_multinom <- as.formula(paste0(target, " ~ ", paste(predictors, collapse="+")))

  dir.create(file.path(output_dir, proj_name_train, "models", "lr-model-multinom", chr), recursive=TRUE)
  dir.create(file.path(output_dir, proj_name_train, "models", "lr-model-summary", chr), recursive=TRUE)
  dir.create(file.path(output_dir, proj_name_train, "models", "lr-model-runtime", chr), recursive=TRUE)

  model_save = file.path(output_dir, proj_name_train, "models", "lr-model-multinom", chr, sprintf("lr_model_k%s_global_val_%s_s%s.rds", k, ct, sample_size))
  summ_save = file.path(output_dir, proj_name_train, "models", "lr-model-summary", chr, sprintf("model_coef_k%s_%s_s%s.rds", k, ct, sample_size))
  runtime_save = file.path(output_dir, proj_name_train, "models", "lr-model-runtime", chr, sprintf("proc_time_k%s_%s_s%s.txt", k, ct, sample_size))
  
  train_1hot_df[,target] <- as.factor(train_1hot_df[,target])
  y_names <- levels(train_1hot_df[target])  # Ensure y is a factor
  levels(train_1hot_df[,target]) <- paste0("X", unlist(lapply(levels(train_1hot_df[,target]), function(x){strsplit(x, split="_")[[1]][1]})))

  y <-  train_1hot_df[,target]
  # Convert data to matrix format (glmnet requires numeric input)
  X <- model.matrix(formula_multinom, data = train_1hot_df)[, -1]  # Remove intercept column

  y_table <- table(y)
  valid_classes <- names(y_table[y_table > 1])  # Keep classes with at least 2 samples
  filtered_indices <- y %in% valid_classes
  X <- X[filtered_indices, ]
  y <- y[filtered_indices]
  y <- droplevels(y)
  print(table(y))


  if((!file.exists(model_save)) | opt$overwrite){
    ptm <- proc.time()
    multinom_model <- glmnet(X, y,family = "multinomial",alpha = 1, lambda=c(0.1,0.01,0.001,0.0001), trace.it=1, standardize=FALSE, type.multinomial="grouped")
    
    ## Get fitted values and save if needed
    if(opt$save_train == TRUE){
      dir.create(file.path(output_dir, proj_name_train, "models", "model-pred-in-train", chr), recursive=TRUE)
      train_pred_save <- file.path(output_dir, proj_name_train, "models", "model-pred-in-train", chr, sprintf("model_pred_train_k%s_%s_s%s.rds", k, ct, sample_size))
      pred_labels <- predict(multinom_model, newx = X, s=0.001, type = "class")
      pred_labels <- factor(pred_labels, levels=paste0("X",seq(1,nstates)))
      levels(pred_labels) <- states_list
      pred_labels <- as.character(pred_labels)

      saveRDS(pred_labels, train_pred_save)

      truth <- train_1hot_df[,target]
      acc <- sum(pred_labels == truth)/length(truth)

      dir.create(file.path(output_dir, proj_name_train, "models", "acc-save", chr), recursive=TRUE)
      capture.output(print(paste0("Accuracy in training celltype ", ct, " is ", acc, "for Lambda = 0.001")), type="output", file=file.path(output_dir, proj_name_train, "models", "acc-save", chr, sprintf("accuracy_k%s_%s_s%s.txt", k, ct, sample_size)))

    }

    ## Save model coef
    model_summ <- coef(multinom_model)
    saveRDS(model_summ, summ_save)

    saveRDS(multinom_model, model_save)
    capture.output(print(paste0("glmnet: ", proc.time()-ptm)), file=runtime_save)
  }
}

print("Models generated!!")


