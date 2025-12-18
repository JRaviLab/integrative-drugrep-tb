# functions to generate intermediate data objects for rank aggregation (for results from both indiv and aggr signatures)
# last modified: 07/10/24
# Kewalin Samart

library(tidyverse)
library(readr)
library(here)
library(dplyr)

source(here("scripts/03_summarize_drugs_methodswise_functions.R"))

rank_drugs_by_connectivity <- function(technologies = c("microarray", "RNAseq"),
                                       scores = c("CMAP", "WCS", "NCS", "Tau", "Cor_spearman", "Cor_pearson"),
                                       run_info_dir = "data/signatures",
                                       results_dir = "results") {
  require(tidyverse)
  require(here)

  score_method_map <- list(
    CMAP = "CMAP",
    WCS = "LINCS",
    NCS = "LINCS",
    Tau = "LINCS",
    Cor_spearman = "Cor",
    Cor_pearson = "Cor"
  )

  for (technology in technologies) {
    message("Processing technology: ", technology)
    data_to_run <- read_tsv(here(run_info_dir, paste0(technology, "_TB_signature_run_info.tsv")))
    data_to_run <- data_to_run[data_to_run$signature == 1, ]

    for (score in scores) {
      score_method <- score_method_map[[score]]
      drug_res_path <- here(results_dir, technology, score_method)

      combined_drug_df <- get_drug_results(data_to_run, drug_res_path, score_method, score)

      score_column <- switch(score,
                             "CMAP" = "scaled_score",
                             "WCS" = "WTCS",
                             "NCS" = "NCS",
                             "Tau" = "Tau",
                             "Cor_spearman" = "cor_score",
                             "Cor_pearson" = "cor_score")

      min_score_df <- combined_drug_df %>%
        dplyr::select(unique_pert = pert, cell, score = !!sym(score_column)) %>%
        group_by(unique_pert, cell) %>%
        summarise(median_score = median(score, na.rm = TRUE), .groups = "drop") %>%
        group_by(unique_pert) %>%
        filter(median_score == min(median_score, na.rm = TRUE)) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        dplyr::rename(min_score = median_score)

      score_matrix <- get_DrugDis_ScoreMatrix(combined_drug_df, score = score, stats = "min")

      ranked_matrix <- score_matrix
      rownames(ranked_matrix) <- ranked_matrix$unique_pert
      ranked_matrix <- ranked_matrix[, -1]

      median_min_scores <- apply(ranked_matrix, 1, median, na.rm = TRUE)
      median_min_scores <- sort(median_min_scores)
      ranked_drugs <- names(median_min_scores)

      min_score_df <- min_score_df[match(ranked_drugs, min_score_df$unique_pert), ]

      median_score_df <- data.frame(
        unique_pert = ranked_drugs,
        median_min_score = as.numeric(median_min_scores),
        stringsAsFactors = FALSE
      )

      output_obj <- list(
        ranked_pert = ranked_drugs,
        min_score_df = min_score_df,
        median_score_df = median_score_df
      )

      output_dir <- here(results_dir, technology, "04_rank_aggregation")
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

      output_file <- file.path(output_dir, paste0(score, "_indivSig_TB_", technology, ".rds"))
      saveRDS(output_obj, output_file)
    }
  }
}

rank_individual_signature_drugs <- function(technologies = c("microarray", "RNAseq"),
                                            scores = c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")) {
  require(here)
  for (technology in technologies) {
    dirname <- here("results", technology, "04_rank_aggregation")
    i <- 1
    for (score in scores) {
      drugs_score_mean_med_added <- readRDS(file.path(dirname, paste0(score, "_indivSig_TB_", technology, ".rds")))
      score_ranked_drugs <- as.data.frame(drugs_score_mean_med_added$ranked_pert)
      colnames(score_ranked_drugs)[1] <- score

      if (i > 1) {
        score_ranked_drugs_df <- cbind(score_ranked_drugs_df, score_ranked_drugs)
      } else {
        score_ranked_drugs_df <- score_ranked_drugs
      }
      i <- i + 1
    }
    saveRDS(score_ranked_drugs_df, file.path(dirname, paste0("full_ranked_drugs_indivSig_TB_", technology, ".rds")))
  }
}

generate_topdrug_lists_from_percentile <- function(technologies = c("microarray", "RNAseq"),
                                                   scores = c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")) {
  require(here)
  for (technology in technologies) {
    dirname <- here("results", technology, "04_rank_aggregation")
    score_ranked_drugs_list <- list()
    top_pct_drugs_list <- list()

    top_drugs_tech_df <- read.delim(
      here("results", technology, "03_methodwise", paste0(technology, "_indiv_top_drugs.tsv")),
      sep = "\t"
    )

    for (score in scores) {
      drugs_score_mean_med_added <- readRDS(here(dirname, paste0(score, "_indivSig_TB_", technology, ".rds")))

      topdrugs_score_mean_med_added <- list(
        ranked_pert = drugs_score_mean_med_added$ranked_pert[
          drugs_score_mean_med_added$ranked_pert %in% top_drugs_tech_df$significant_drug
        ],
        min_score_df = drugs_score_mean_med_added$min_score_df[
          drugs_score_mean_med_added$min_score_df$unique_pert %in% top_drugs_tech_df$significant_drug, ]
      )

      score_ranked_drugs <- topdrugs_score_mean_med_added$ranked_pert
      top_pct_drugs_list[[score]] <- top_drugs_tech_df$ranked_pert
      score_ranked_drugs_list[[score]] <- score_ranked_drugs
    }

    saveRDS(score_ranked_drugs_list, file.path(dirname, paste0("ranked_topdrugs_indivSig_TB_", technology, "_list.rds")))
  }
}

rank_aggregated_signature_drugs <- function(technologies = c("microarray", "RNAseq"),
                                            score_order = c("CMAP", "WCS", "NCS", "Tau", "Cor_spearman", "Cor_pearson"),
                                            values = "neg",
                                            tail = "left",
                                            score_percentile = 0.9,
                                            base_dir = "results") {
  require(here)
  require(readr)
  require(dplyr)

  suffix <- paste0(values, "_", tail, "_", score_percentile)

  for (technology in technologies) {
    message("Processing: ", technology)

    methodwise_dir <- here(base_dir, technology, "03_methodwise")
    top_drug_file <- file.path(methodwise_dir, paste0(technology, "_aggr_top_drugs.tsv"))

    if (!file.exists(top_drug_file)) {
      warning("Top drug file missing: ", top_drug_file)
      next
    }

    top_drugs <- read_tsv(top_drug_file, show_col_types = FALSE)$significant_drug
    drug_lists <- list()

    for (score in score_order) {
      score_file <- file.path(methodwise_dir, paste0(score, "_aggrSig_ranked_", suffix, ".tsv"))

      if (!file.exists(score_file)) {
        warning("Score file not found: ", score_file)
        next
      }

      df <- read_tsv(score_file, show_col_types = FALSE)

      if (!"pert" %in% colnames(df)) {
        warning("Missing 'pert' column in: ", score_file)
        next
      }

      # Filter ranked drugs by those in top_drugs
      ranked <- df$pert[df$pert %in% top_drugs]
      drug_lists[[score]] <- ranked
    }

    # Get unique top drugs actually ranked in at least one method
    filtered_top_drugs <- unique(unlist(drug_lists))

    # Initialize rank matrix
    rank_matrix <- matrix(NA, nrow = length(filtered_top_drugs), ncol = length(score_order))
    rownames(rank_matrix) <- filtered_top_drugs
    colnames(rank_matrix) <- score_order

    for (score in score_order) {
      if (!is.null(drug_lists[[score]])) {
        rank_matrix[, score] <- match(filtered_top_drugs, drug_lists[[score]])
      }
    }

    score_ranked_drugs_df <- as.data.frame(rank_matrix)
    score_ranked_drugs_df$pert <- rownames(score_ranked_drugs_df)
    score_ranked_drugs_df <- score_ranked_drugs_df[, c("pert", score_order)]

    # Convert to list of ordered drug vectors
    score_ranked_drugs_list <- lapply(score_order, function(score) {
      score_ranked_drugs_df[order(score_ranked_drugs_df[[score]], na.last = NA), "pert"]
    })
    names(score_ranked_drugs_list) <- score_order

    # Save output
    save_dir <- file.path(here(base_dir, technology, "04_rank_aggregation"))
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

    saveRDS(score_ranked_drugs_df,
            file.path(save_dir, paste0("full_ranked_drugs_aggrSig_TB_", technology, ".rds")))
    print(paste0("Succesfully saved full_ranked_drugs_aggrSig_TB_", technology, ".rds", " at ", save_dir))

    saveRDS(score_ranked_drugs_list,
            file.path(save_dir, paste0("ranked_aggrSig_drugs_scores_", suffix, "_", technology, "_list.rds")))
    print(paste0("Succesfully saved ranked_aggrSig_drugs_scores_", suffix, "_", technology, "_list.rds", " at ", save_dir))
  }
}
