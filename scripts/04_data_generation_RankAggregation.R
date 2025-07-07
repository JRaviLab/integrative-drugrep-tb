# script to generate intermediate data objects for rank aggregation (for results from both indiv and aggr signatures)
# created date: 01/02/24
# last modified: 07/07/24
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
data_to_run <- read_tsv(here("data/v2/signatures/microarray_TB_signature_run_info.tsv"))  # Adjust path if needed
data_to_run <- data_to_run[data_to_run$signature == 1,]

# Main loop
for (technology in technologies) {
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

    # Get min score + cell type per drug
    min_score_df <- combined_drug_df %>%
      select(unique_pert = pert, cell, score = !!sym(score_column)) %>%
      group_by(unique_pert) %>%
      filter(score == min(score, na.rm = TRUE)) %>%
      slice(1) %>%
      ungroup() %>%
      rename(min_score = score)

    # Get score matrix (rows: drugs, cols: signatures)
    score_matrix <- get_DrugDis_ScoreMatrix(combined_drug_df, score = score, stats = "min")

    # Rank each column, then aggregate ranks across signatures
    ranked_matrix <- score_matrix
    rownames(ranked_matrix) <- ranked_matrix$unique_pert
    ranked_matrix <- ranked_matrix[, -1]

    rank_df <- apply(ranked_matrix, 2, function(x) rank(x, na.last = "keep"))
    mean_ranks <- rowMeans(rank_df, na.rm = TRUE)
    ranked_drugs <- names(sort(mean_ranks))
    # Reorder min_score_df to match ranked_drugs
    min_score_df <- min_score_df[match(ranked_drugs, min_score_df$unique_pert), ]

    # Save as RDS
    output_obj <- list(
      ranked_pert = ranked_drugs,
      min_score_df = min_score_df  # includes unique_pert, min_score, cell
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
for(data_technology in c("microarray","RNAseq")){
  dirname <- paste0("results",data_technology,"/")
  for(score in scores){
    top_drugs_tech_df <- read.delim(paste0("./results/uniformly_processed/",data_technology,"/",score,"_top_0.9pct_0.5reversed_drugs.tsv"),sep = "\t")
    drugs_score_mean_med_added <- readRDS(paste0(dirname,score,"_indivSig_TB_",data_technology,"_matrix.rds")) # overwrite the current matrix with same data and mean/median added
    topdrugs_score_mean_med_added <- drugs_score_mean_med_added[drugs_score_mean_med_added$unique_pert %in% top_drugs_tech_df$unique_pert,]
    score_ranked_drugs <- topdrugs_score_mean_med_added$unique_pert
    top_pct_drugs_list[[score]] <- top_drugs_tech_df$unique_pert
    score_ranked_drugs_list[[score]] <- score_ranked_drugs
  }
  saveRDS(score_ranked_drugs_list, paste0(dirname,"04_rank_aggregation/ranked_0.9pct0.5rv_topdrugs_indivSig_TB_",data_technology,"_list.rds"))
  saveRDS(top_pct_drugs_list, paste0(dirname,"04_rank_aggregation/unranked_0.9pct0.5rv_topdrugs_indivSig_TB_",data_technology,"_list.rds"))
}

#-----------------
# Aggregated signatures
# generate lists of partial ranked top drugs
for(data_technology in c("microarray","RNAseq")){
  dirname <- paste0("./results/uniformly_processed/",data_technology,"/")
  drugs_score <- read.csv(paste0(dirname,"uni_",data_technology,"_aggrSig_drugs_scores_neg_left_0.9.tsv"),sep="\t")
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

