# script to generate intermediate data objects for rank aggregation (for results from both indiv and aggr signatures)
# created date: 01/02/24
# last modified: 07/08/24
# Kewalin Samart

library(tidyverse)
library(readr)
library(here)
library(dplyr)

source(here("scripts/03_summarize_drugs_methodswise_functions.R"))

# Define variables
technologies <- c("microarray", "RNAseq")
scores <- c("CMAP", "WCS", "NCS", "Tau", "Cor_spearman", "Cor_pearson")
score_method_map <- list(
  CMAP = "CMAP",
  WCS = "LINCS",
  NCS = "LINCS",
  Tau = "LINCS",
  Cor_spearman = "Cor",
  Cor_pearson = "Cor"
)

# Load metadata table that contains SIGNATURE_NAME


# Main loop
for (technology in technologies) {
  print(technology)
  data_to_run <- read_tsv(here(paste0("data/v2/signatures/",technology,"_TB_signature_run_info.tsv")))  # Adjust path if needed
  data_to_run <- data_to_run[data_to_run$signature == 1,]
  for (score in scores) {
    score_method <- score_method_map[[score]]
    drug_res_path <- here(paste0("results/", technology,"/",score_method))

    # Get combined drug results for all signatures
    combined_drug_df <- get_drug_results(data_to_run, drug_res_path, score_method, score)

    # Identify correct score column
    score_column <- switch(score,
                           "CMAP" = "scaled_score",
                           "WCS" = "WTCS",
                           "NCS" = "NCS",
                           "Tau" = "Tau",
                           "Cor_spearman" = "cor_score",
                           "Cor_pearson" = "cor_score"
    )
    # the cell line where it shows the lowest median score across all disease signatures
    min_score_df <- combined_drug_df %>%
      select(unique_pert = pert, cell, score = !!sym(score_column)) %>%
      group_by(unique_pert, cell) %>%
      summarise(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
      group_by(unique_pert) %>%
      filter(median_score == min(median_score, na.rm = TRUE)) %>%
      slice(1) %>%
      ungroup() %>%
      rename(min_score = median_score)

    # Get score matrix (rows: drugs, cols: signatures)
    score_matrix <- get_DrugDis_ScoreMatrix(combined_drug_df, score = score, stats = "min")

    # Rank each column, then aggregate ranks across signatures
    ranked_matrix <- score_matrix
    rownames(ranked_matrix) <- ranked_matrix$unique_pert
    ranked_matrix <- ranked_matrix[, -1]

    # Compute median of minimum scores across all signatures (rows = drugs)
    median_min_scores <- apply(ranked_matrix, 1, median, na.rm = TRUE)

    # Sort drugs: smaller median min = stronger reversal → rank 1
    median_min_scores <- sort(median_min_scores)  # ascending order
    ranked_drugs <- names(median_min_scores)

    # Reorder min_score_df to match new ranked order
    min_score_df <- min_score_df[match(ranked_drugs, min_score_df$unique_pert), ]

    # the median of minimum scores across all disease signatures, regardless of cell line
    median_score_df <- data.frame(
      unique_pert = ranked_drugs,
      median_min_score = as.numeric(median_min_scores),
      stringsAsFactors = FALSE
    )

    # Add median_min_scores as a data frame for output
    median_score_df <- data.frame(
      unique_pert = names(median_min_scores),
      median_min_score = as.numeric(median_min_scores),
      stringsAsFactors = FALSE
    )

    # Reorder min_score_df to match ranked_drugs
    min_score_df <- min_score_df[match(ranked_drugs, min_score_df$unique_pert), ]

    # Save as RDS
    output_obj <- list(
      ranked_pert = ranked_drugs,
      min_score_df = min_score_df,  # includes best score per drug
      median_score_df = median_score_df  # NEW: includes median score per drug
    )


    output_dir <- here("results", technology, "04_rank_aggregation")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    output_file <- file.path(output_dir, paste0(score, "_indivSig_TB_", technology, ".rds"))
    saveRDS(output_obj, output_file)
  }
}

#----------------------------------------

# Individual signatures
# generate full ranked drugs matrix -- individual signatures
scores <- c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")
# compute mean and median of each drug
for(technology in c("microarray","RNAseq")){
  dirname <- here(paste0("results/",technology,"/04_rank_aggregation/"))
  i = 1
  for(score in scores){
    drugs_score_mean_med_added <- readRDS(paste0(dirname,score,"_indivSig_TB_",technology,".rds"))
    score_ranked_drugs <- as.data.frame(drugs_score_mean_med_added$ranked_pert)
    colnames(score_ranked_drugs)[1] <- score
    if(i > 1){
      score_ranked_drugs_df <- cbind(score_ranked_drugs_df,score_ranked_drugs)
    }else{
      score_ranked_drugs_df <- score_ranked_drugs
    }
    i = i + 1
  }
  saveRDS(score_ranked_drugs_df, paste0(dirname,"/full_ranked_drugs_indivSig_TB_",technology,".rds"))
}

#--------------------------------------------------------------------

# generate lists of partial ranked top drugs and unranked top drugs from pct approach
top_pct_drugs_list <- list()
score_ranked_drugs_list <- list()
scores <- c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")
for(technology in c("microarray","RNAseq")){
  dirname <- paste0("results/",technology,"/04_rank_aggregation/")
  for(score in scores){
    top_drugs_tech_df <- read.delim(here(paste0("results/",technology,"/03_methodwise/",technology,"_indiv_top_drugs.tsv")),sep = "\t")
    drugs_score_mean_med_added <- readRDS(here(paste0(dirname,score,"_indivSig_TB_",technology,".rds"))) # overwrite the current matrix with same data and mean/median added
    #topdrugs_score_mean_med_added <- drugs_score_mean_med_added[drugs_score_mean_med_added$ranked_pert %in% top_drugs_tech_df$significant_drug,]
    topdrugs_score_mean_med_added <- list(
      ranked_pert = drugs_score_mean_med_added$ranked_pert[drugs_score_mean_med_added$ranked_pert %in% top_drugs_tech_df$significant_drug],
      min_score_df = drugs_score_mean_med_added$min_score_df[
        drugs_score_mean_med_added$min_score_df$unique_pert %in% top_drugs_tech_df$significant_drug, ]
    )
    score_ranked_drugs <- topdrugs_score_mean_med_added$ranked_pert
    top_pct_drugs_list[[score]] <- top_drugs_tech_df$ranked_pert
    score_ranked_drugs_list[[score]] <- score_ranked_drugs
  }
  saveRDS(score_ranked_drugs_list, paste0(dirname,"ranked_topdrugs_indivSig_TB_",technology,"_list.rds"))
}

#-----------------
# Aggregated signatures
# generate lists of partial ranked top drugs
for(technology in c("microarray","RNAseq")){
  dirname <- paste0("results/uniformly_processed/",technology,"/")
  drugs_score <- read.csv(here(paste0(dirname, technology,"_aggrSig_drugs_scores_neg_left_0.9.tsv")),sep="\t")
  for(i in 3:dim(drugs_score)[2]){
    score_ranked_drugs = drugs_score[order(drugs_score[i], decreasing = FALSE), ]$pert
    score_ranked_drugs <- as.data.frame(score_ranked_drugs)
    if(i > 3){
      score_ranked_drugs_df <- cbind(score_ranked_drugs_df,score_ranked_drugs)
    }else{
      score_ranked_drugs_df <- score_ranked_drugs
    }
  }
  colnames(score_ranked_drugs_df) <- c("CMAP","Tau","NCS","WCS","Cor_spearman","Cor_pearson")
  score_ranked_drugs_df <- score_ranked_drugs_df[c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")]
  score_ranked_drugs_list <- lapply(names(score_ranked_drugs_df), function(col) score_ranked_drugs_df[[col]])
  names(score_ranked_drugs_list) <- colnames(score_ranked_drugs_df)
  saveRDS(score_ranked_drugs_df, paste0(dirname,"04_rank_aggregation/full_ranked_drugs_aggrSig_TB_",data_technology,"_matrix.rds"))
  saveRDS(score_ranked_drugs_list, paste0(dirname,"04_rank_aggregation/ranked_aggrSig_drugs_scores_neg_left_0.9_",data_technology,"_list.rds"))
}

