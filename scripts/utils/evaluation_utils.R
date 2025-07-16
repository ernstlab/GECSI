### GENERAL FUNCTIONS 
# Plot the precision-recall curve
pr_curve <- function(interp_pr_data, pr_data, save_dir){
    p1 <- ggplot(interp_pr_data, aes(x = recall, y = precision)) +
    geom_step(direction="vh") +
    geom_point(data = pr_data, aes(x = recall, y = precision), color = "red") +
    labs(title = "Precision-Recall Curve",
        x = "Recall",
        y = "Precision") +
    theme_minimal()

    ggsave(file=save_dir, p1, width=6, height=6)
}

### BASIC HELPER FUNCTIONS ###
set_colors <- function(df, col, colors="default"){
    if(colors == "default"){
        n <- length(unique(df[,col]))
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        set.seed(10)
        cols=sample(col_vector, n)
    } else {
        cols=colors
    }
    return(cols)
}

calculate_precision_recall <- function(predictions, labels) {
  # Sort predictions and corresponding labels
  sorted_indices <- order(predictions, decreasing = TRUE)
  sorted_predictions <- predictions[sorted_indices]
  sorted_labels <- labels[sorted_indices]
  
  # Calculate true positives and false positives using cumsum
  true_positives <- cumsum(sorted_labels)
  false_positives <- cumsum(!sorted_labels)
  
  # Number of positives in the labels
  n_positives <- sum(labels == 1)
  
  # Calculate recall and precision vectors
  recall <- true_positives / n_positives
  precision <- true_positives / (true_positives + false_positives)
  
  # Remove duplicates in recall
  unique_recall <- !duplicated(recall)
  
  # Only keep unique precision-recall pairs
  precision <- precision[unique_recall]
  recall <- recall[unique_recall]
  thresholds <- sorted_predictions[unique_recall]
  
  return(data.frame(threshold = thresholds, precision = precision, recall = recall))
}

# calculate_precision_recall <- function(predictions, labels) {
#   thresholds <- unique(predictions)
#   thresholds <- sort(thresholds, decreasing = TRUE)
  
#   precision <- numeric(length(thresholds))
#   recall <- numeric(length(thresholds))
  
#   for (i in seq_along(thresholds)) {
#     threshold <- thresholds[i]
#     predicted_positive <- predictions >= threshold
    
#     true_positive <- sum(predicted_positive & labels == 1)
#     false_positive <- sum(predicted_positive & labels == 0)
#     false_negative <- sum(!predicted_positive & labels == 1)
    
#     precision[i] <- ifelse((true_positive + false_positive) == 0, 1, true_positive / (true_positive + false_positive))
#     recall[i] <- true_positive / (true_positive + false_negative)
#   }
  
#   return(data.frame(threshold = thresholds, precision = precision, recall = recall))
# }

interpolate_pr_curve_and_auc <- function(pr_data) {
  interp_points <- data.frame(
    recall = c(0, pr_data$recall, 1),
    precision = c(1, pr_data$precision, 0)
  )
  # Calculate the area under the curve (AUC-PR) using trapezoidal rule
  # auc_pr <- sum(diff(interp_points$recall) * (interp_points$precision[-1] + interp_points$precision[-length(interp_points$precision)]) / 2)
  auc_pr <- sum(interp_points$precision[-1]*(interp_points$recall[-1] - interp_points$recall[-length(interp_points$recall)]))  
  return(list(interpolated_points = interp_points, auc_pr = auc_pr))
}

add_locus_dist_info <- function(df_out, df_src){
    if(nrow(df_out) == nrow(df_src)){
       df_out$locus <- df_src$locus
       df_out$Dist1 <- df_src$Dist1
       return(df_out)
    } else {
        stop("The number of rows of two data frames are not the same!")
    }
}

make_long <- function(df, ct, value){
  df_long <- pivot_longer(df, cols=all_of(ct), names_to="Cell_type", values_to=value)
  df_long <- as.data.frame(df_long)
  return(df_long)
}

read_wig <- function(x, format='wig', genome='mm9') {
  # keepStandardChromosomes()
  suppressMessages(library(rtracklayer))  
  merged_wig <- import.wig(x, format=format, genome=genome)
  merged_wig <- keepSeqlevels(merged_wig, paste0('chr', c(seq(1,19), 'X', 'Y')), pruning.mode="coarse")
  return(merged_wig)
}

extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*")){ warning("'*' ranges were treated as '+'") }
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  # trim(x)
  x
}

reduce_to_tss <- function(x) {
  return( resize(x, 1, fix="start", use.names=TRUE, ignore.strand=FALSE) ) 
}

mround <- function(x,base){ 
  base*round(x/base) 
} 


### UTIL FUNCTIONS FOR GETTING FREQUENCY
most_frequent_state_in_row <- function(row) {
  frequency_table <- table(unname(unlist(row)))
  most_frequent_state <- names(frequency_table[which.max(frequency_table)])
  return(most_frequent_state)
}
