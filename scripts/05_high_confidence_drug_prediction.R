# Getting high confidence drug candidates
# last modified: 12/18/25
# Author: Ling Thang
# Kewalin Samart edited

# import libraries
library(here)
library(dplyr)
library(readr)

# function to read rank-aggregated results from four pipelines with input validation
read_pipeline_results <- function(microarray_indiv_path, microarray_agg_path, rnaseq_indiv_path, rnaseq_agg_path) {
  # helper to check file existence
  check_file <- function(path) {
    if (!file.exists(path)) {
      stop(paste("File does not exist:", path))
    }
  }
  check_file(microarray_indiv_path)
  check_file(microarray_agg_path)
  check_file(rnaseq_indiv_path)
  check_file(rnaseq_agg_path)

  # read files
  microarray_indiv <- read_tsv(microarray_indiv_path, show_col_types = FALSE)
  microarray_agg <- read_tsv(microarray_agg_path, show_col_types = FALSE)
  rnaseq_indiv <- read_tsv(rnaseq_indiv_path, show_col_types = FALSE)
  rnaseq_agg <- read_tsv(rnaseq_agg_path, show_col_types = FALSE)

  # validate required columns
  required_cols <- c("drug_name", "rank_score")
  for (df in list(microarray_indiv, microarray_agg, rnaseq_indiv, rnaseq_agg)) {
    missing <- setdiff(required_cols, colnames(df))
    if (length(missing) > 0) {
      stop(paste("Missing columns:", paste(missing, collapse = ", ")))
    }
  }

  return(list(
    indivSig_TB_microarray = microarray_indiv,
    aggSig_TB_microarray = microarray_agg,
    indivSig_TB_RNAseq = rnaseq_indiv,
    aggSig_TB_RNAseq = rnaseq_agg
  ))
}

# process individual and aggregated pipelines separately
process_pipeline <- function(indiv_microarray, indiv_rnaseq, agg_microarray, agg_rnaseq) {
  # individual signatures pipeline
  indiv_merged_df <- dplyr::full_join(indiv_microarray, indiv_rnaseq, by = "drug_name")
  indiv_results <- indiv_merged_df %>%
    rowwise() %>%
    mutate(mean_rank_score_indiv = mean(c(rank_score.x, rank_score.y), na.rm = TRUE)) %>%
    dplyr::select(drug_name, mean_rank_score_indiv)
  message("Processed individual signatures pipeline.")

  # aggregated signatures pipeline
  agg_merged_df <- dplyr::full_join(agg_microarray, agg_rnaseq, by = "drug_name")
  agg_results <- agg_merged_df %>%
    rowwise() %>%
    mutate(mean_rank_score_agg = mean(c(rank_score.x, rank_score.y), na.rm = TRUE)) %>%
    dplyr::select(drug_name, mean_rank_score_agg)
  message("Processed aggregated signatures pipeline.")

  return(list(indiv_results = indiv_results, agg_results = agg_results))
}

convert_to_zscore <- function(scores_vector) {
  #' @description Converts a numeric vector into standard scores (z-scores).
  #' @param scores_vector A numeric vector of scores.
  #' @returns A numeric vector of z-scores.

  z_scores <- (scores_vector - mean(scores_vector, na.rm = TRUE)) / sd(scores_vector, na.rm = TRUE)
  return(z_scores)
}

# identify high-confidence drugs and calculate final score
get_high_confidence_predictions <- function(indiv_results, agg_results) {
  final_merged_df <- merge(indiv_results, agg_results, by = "drug_name", all = FALSE)

  message(paste("Identified", nrow(final_merged_df), "combined drug candidates."))

  # convert each pipeline's mean rank score to a z-score
  final_merged_df$z_score_indiv <- convert_to_zscore(final_merged_df$mean_rank_score_indiv)
  final_merged_df$z_score_agg <- convert_to_zscore(final_merged_df$mean_rank_score_agg)

  # calculate the combined z-score as per the formula in Figure 1c
  combined_predictions <- final_merged_df %>%
    mutate(
      # formula: (z1 + z2) / sqrt(2)
      combined_z_score = (z_score_indiv + z_score_agg) / sqrt(2)
    ) %>%
    dplyr::select(drug_name, combined_z_score) %>%
    arrange(desc(combined_z_score))

  # calculate the final combined z-score as per the formula in Figure 1c
  high_confidence_predictions <- final_merged_df %>%
    mutate(
      # formula: (z1 + z2) / sqrt(2)
      combined_z_score = (z_score_indiv + z_score_agg) / sqrt(2)
    ) %>%
    dplyr::filter(z_score_indiv > 0, z_score_agg > 0, combined_z_score > 0) %>%
    dplyr::select(drug_name, combined_z_score) %>%
    arrange(desc(combined_z_score))
  message(paste("Identified", nrow(high_confidence_predictions), "high confidence drug candidates."))

  return(list(combined_predictions, high_confidence_predictions))
}

# save and display results
save_and_display_predictions <- function(predictions, output_path) {
  write_tsv(predictions, output_path)
  message(paste("Final high-confidence drug predictions saved to:", output_path))
  print("Top 20 High-Confidence Drug Predictions:")
  print(head(predictions, 20))
}
