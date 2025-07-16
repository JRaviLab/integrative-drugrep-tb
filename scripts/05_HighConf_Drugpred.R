# utility functions for #drugrep_tb project
# created date: 07/16/25
# last modified: 07/16/25
# Ling Thang (edited by Gemini)

# Implements Figure 1c from Integrative transcriptome-based drug repurposing in tuberculosis
# Key implementation (page 10) :
# "(c) Drug prioritization pipeline applied to microarray and RNA-seq signatures, including both individual and aggregated
# versions. Drugs reversing ≥50% of the disease signature with scores in the top 10% most negative values were selected,
# and only those supported by ≥2 of 3 metric subgroups (CMAP 1.0, LINCS, correlation-based) were retained. Final drug rankings
# were generated using Bayesian Latent Variable Approach for partial rank aggregation; “BiG” method [44], and high-confidence predictions
# were obtained by intersecting results across signature types."

# 1. Setup

# import libraries
library(here)
library(dplyr)
library(readr)

# Define paths to the rank-aggregated results from the four pipelines
# MICROARRAY
indivSig_TB_microarray_path <- here("results/microarray/04_rank_aggregation/RAresult_indivSig_TB_microarray.tsv")
aggSig_TB_microarray_path <- here("results/microarray/04_rank_aggregation/RAresult_neg_left_0.9_aggrSig_TB_microarray.tsv")

# RNA-SEQ
indivSig_TB_RNAseq_path <- here("results/RNAseq/04_rank_aggregation/RAresult_indivSig_TB_RNAseq.tsv")
aggSig_TB_RNAseq_path <- here("results/RNAseq/04_rank_aggregation/RAresult_neg_left_0.9_aggrSig_TB_RNAseq.tsv")

# Read data into dataframes 'drug_name' and 'rank_score' !EXPECTED!
indivSig_TB_microarray <- read_tsv(indivSig_TB_microarray_path, show_col_types = FALSE)
aggSig_TB_microarray <- read_tsv(aggSig_TB_microarray_path, show_col_types = FALSE)

indivSig_TB_RNAseq <- read_tsv(indivSig_TB_RNAseq_path, show_col_types = FALSE)
aggSig_TB_RNAseq <- read_tsv(aggSig_TB_RNAseq_path, show_col_types = FALSE)

# 2. Utility Function

convert_to_zscore <- function(scores_vector) {
    #' @description Converts a numeric vector into standard scores (z-scores).
    #' @param scores_vector A numeric vector of scores.
    #' @returns A numeric vector of z-scores.

    z_scores <- (scores_vector - mean(scores_vector, na.rm = TRUE)) / sd(scores_vector, na.rm = TRUE)
    return(z_scores)
}

# 3. Process Individual and Aggregated Pipelines Separately

# --- Individual Signatures Pipeline ---
# Use a full join to merge the two dataframes.
# dplyr automatically handles the identical 'rank_score' column by adding suffixes (.x, .y)
indiv_merged_df <- dplyr::full_join(indivSig_TB_microarray, indivSig_TB_RNAseq, by = "drug_name")

# Calculate the mean rank score for the individual pipeline using the new column names
indiv_results <- indiv_merged_df %>%
    rowwise() %>%
    mutate(mean_rank_score_indiv = mean(c(rank_score.x, rank_score.y), na.rm = TRUE)) %>%
    dplyr::select(drug_name, mean_rank_score_indiv)

message("Processed individual signatures pipeline.")


# --- Aggregated Signatures Pipeline ---
# Use a full join to merge the two dataframes.
agg_merged_df <- dplyr::full_join(aggSig_TB_microarray, aggSig_TB_RNAseq, by = "drug_name")

# Calculate the mean rank score for the aggregated pipeline
agg_results <- agg_merged_df %>%
    rowwise() %>%
    mutate(mean_rank_score_agg = mean(c(rank_score.x, rank_score.y), na.rm = TRUE)) %>%
    dplyr::select(drug_name, mean_rank_score_agg)

message("Processed aggregated signatures pipeline.")

# 4. Identify High-Confidence Drugs and Calculate Final Score

# Merge the results from both pipelines (inner join automatically finds high-confidence drugs)
final_merged_df <- merge(indiv_results, agg_results, by = "drug_name", all = FALSE)

message(paste("Identified", nrow(final_merged_df), "high-confidence drug candidates."))

# Convert each pipeline's mean rank score to a z-score
final_merged_df$z_score_indiv <- convert_to_zscore(final_merged_df$mean_rank_score_indiv)
final_merged_df$z_score_agg <- convert_to_zscore(final_merged_df$mean_rank_score_agg)

# Calculate the final combined z-score as per the formula in Figure 1c
high_confidence_predictions <- final_merged_df %>%
    mutate(
        # Formula: (z1 + z2) / sqrt(2)
        combined_z_score = (z_score_indiv + z_score_agg) / sqrt(2)
    ) %>%
    dplyr::select(drug_name, combined_z_score) %>%
    arrange(desc(combined_z_score))

# 5. Save and Display Results

output_path <- here("results/final_high_confidence_drug_predictions.tsv")
write_tsv(high_confidence_predictions, output_path)

message(paste("Final high-confidence drug predictions saved to:", output_path))

# Display the top 20 predicted drugs
print("Top 20 High-Confidence Drug Predictions:")
print(head(high_confidence_predictions, 20))
