# script to perform rank aggregation on the ranked drug results
# Used approach: BiG (A Bayesian latent variable approach i.e., Bayesian Aggregation in Genomic applications)
# **for partial and top ranked lists; incorporating the effect of clustering**
# ref1: https://doi.org/10.1002%2Fsim.7920
# ref2: https://github.com/baillielab/comparison_of_RA_methods/blob/main/algorithms/useBiG.R
# created date: 01/02/24
# last modified: 07/08/25
# Kewalin Samart

library(tidyverse)
library(readr)
library(here)

source(here("scripts/04_RankAggregation_utilities.R"))
source(here("scripts/04_BiG_RankAggregation_functions.R"))

run_BiG_diffuse_indivSig <- function(data_technologies = c("microarray", "RNAseq"),
                                                 n_p1_params = c(4, 4),
                                                 n_iter = 2000,
                                                 burnin = 1000,
                                                 prior = "IG") {
  require(here)
  require(readr)

  for (i in seq_along(data_technologies)) {
    data_technology <- data_technologies[i]
    n_p1 <- n_p1_params[i]

    dirname <- here("results", data_technology, "04_rank_aggregation")

    # load ranked top drug list
    ranked_drug_list_path <- file.path(dirname, paste0("ranked_topdrugs_indivSig_TB_", data_technology, "_list.rds"))
    ranked_drug_list <- readRDS(ranked_drug_list_path)

    # convert to rank matrix
    rank_matrix <- matrix_transfer(ranked_drug_list, full = FALSE)

    # determine number of top-ranked drugs per score
    NTlength <- sapply(rank_matrix, length)

    # run BiG diffuse aggregation
    result <- BiG_diffuse(r = rank_matrix, n_T = NTlength, n_p1 = n_p1,
                          M = n_iter, burnin = burnin, prior = prior)

    # extract ranked drug entities
    entities <- rownames(rank_matrix)
    rankedEntities <- entities[order(result, decreasing = TRUE)]

    # create and save data frame of ranked scores
    RA_df <- data.frame(drug_name = entities, rank_score = result)
    RA_df <- RA_df[order(-RA_df$rank_score), ]

    # save outputs
    write_tsv(RA_df, file.path(dirname, paste0("RAresult_indivSig_TB_", data_technology, ".tsv")))
    saveRDS(rankedEntities, file.path(dirname, paste0("RAresult_indivSig_TB_", data_technology, ".rds")))
  }
}

run_BiG_diffuse_aggrSig <- function(technologies = c("microarray", "RNAseq"),
                                       n_p1_params = c(4, 4),
                                       base_dir = "results",
                                       values = "neg",
                                       tail = "left",
                                       score_percentile = 0.9,
                                       M = 2000,
                                       burnin = 1000,
                                       prior = "IG") {
  require(readr)

  suffix <- paste0(values, "_", tail, "_", score_percentile)

  for (i in seq_along(technologies)) {
    technology <- technologies[i]
    n_p1 <- n_p1_params[i]

    message("Running BiG rank aggregation for: ", technology)

    dirname <- file.path(base_dir, technology)
    rank_file <- here(dirname, "04_rank_aggregation",
                           paste0("ranked_aggrSig_drugs_scores_", suffix, "_", technology, "_list.rds"))

    if (!file.exists(here(rank_file))) {
      warning("Ranked drug list missing for: ", technology)
      next
    }

    # Read ranked drug list and transform to rank matrix
    ranked_drug_list <- readRDS(rank_file)
    rank_matrix <- matrix_transfer(ranked_drug_list, full = TRUE)

    # Adjust ranks using 'min' tie method
    rank_matrix_adjusted <- as.data.frame(lapply(as.data.frame(rank_matrix), function(x) rank(x, ties.method = "min")))
    rownames(rank_matrix_adjusted) <- rownames(rank_matrix)
    rank_matrix_adjusted <- as.matrix(rank_matrix_adjusted)

    # Calculate NTlength (optional for BiG, might be overridden)
    NTlength <- sapply(seq_len(ncol(rank_matrix_adjusted)), function(i) length(rank_matrix_adjusted[, i]))

    # Run BiG Diffuse aggregation
    result <- BiG_diffuse(r = rank_matrix_adjusted,
                          n_T = NTlength,
                          n_p1 = n_p1,
                          M = M,
                          burnin = burnin,
                          prior = prior)

    # Rank drugs
    entities <- rownames(rank_matrix_adjusted)
    rankedEntities <- entities[order(result, decreasing = TRUE)]
    RA_result_df <- data.frame(drug_name = entities, rank_score = result) %>%
      arrange(desc(rank_score))

    # Save RA outputs
    save_dir <- here(dirname, "04_rank_aggregation")
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

    write_tsv(RA_result_df,
              file.path(save_dir, paste0("RAresult_", suffix, "_aggrSig_TB_", technology, ".tsv")))
    saveRDS(rankedEntities,
            file.path(save_dir, paste0("RAresult_", suffix, "_aggrSig_TB_", technology, ".rds")))
  }
}
