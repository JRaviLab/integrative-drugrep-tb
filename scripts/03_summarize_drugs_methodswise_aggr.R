# functions to get drug result summary and top drugs from an aggregated signature
# last modified: 07/10/25
# Kewalin Samart

library(readr)
library(here)

source(here("scripts/03_summarize_drugs_methodswise_functions.R"))

extract_aggrSig_significant_drugs <- function(technologies = c("microarray", "RNAseq"),
                                                      scores = c("CMAP", "WCS", "NCS", "Tau", "Cor_spearman", "Cor_pearson"),
                                                      values = "neg",
                                                      tail = "left",
                                                      score_percentile = 0.9,
                                                      base_dir = "results") {
  require(dplyr)
  require(readr)
  require(here)

  score_method_map <- list(
    CMAP = "CMAP",
    WCS = "LINCS",
    NCS = "LINCS",
    Tau = "LINCS",
    Cor_spearman = "Cor",
    Cor_pearson = "Cor"
  )

  signi_drugs_score_df_all <- list()
  signi_drugs_score_list_all <- list()

  for (technology in technologies) {
    for (score in scores) {
      score_method <- score_method_map[[score]]

      file_path <- if (score %in% c("WCS", "NCS", "Tau")) {
        here(base_dir, technology, score_method, "LINCS_aggregated_signature.tsv")
      } else {
        here(base_dir, technology, score_method, paste0(score, "_aggregated_signature.tsv"))
      }

      if (!file.exists(file_path)) {
        warning("Skipping missing file: ", file_path)
        next
      }

      drug_df <- read.delim(file_path, sep = "\t")

      score_col <- switch(score,
                          Tau = "Tau",
                          NCS = "NCS",
                          WCS = "WTCS",
                          CMAP = "scaled_score",
                          Cor_spearman = "cor_score",
                          Cor_pearson = "cor_score"
      )

      score_vec <- drug_df[[score_col]]
      threshold <- determine_threshold(score_vec, values = values, tail = tail, score_percentile = score_percentile)

      if (values == "neg") {
        signi_drug_df <- drug_df[drug_df[[score_col]] < threshold, ]
      } else if (values == "pos") {
        signi_drug_df <- drug_df[drug_df[[score_col]] > threshold, ]
      } else {
        signi_drug_df <- drug_df[abs(drug_df[[score_col]]) < abs(threshold), ]
      }

      signi_drug_df <- signi_drug_df[, c("pert", "cell", score_col)]

      if (score %in% c("Cor_spearman", "Cor_pearson")) {
        colnames(signi_drug_df)[3] <- ifelse(score == "Cor_spearman", "spearman", "pearson")
        score_col <- colnames(signi_drug_df)[3]
      }

      signi_drug_df <- signi_drug_df %>% arrange(!!sym(score_col))

      key <- paste(technology, score, sep = "_")
      signi_drugs_score_df_all[[key]] <- signi_drug_df
      signi_drugs_score_list_all[[key]] <- signi_drug_df$pert

      # Save output
      output_dir <- here(base_dir, technology, "03_methodwise")
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

      out_file <- file.path(output_dir, paste0(score, "_aggrSig_ranked_", values, "_", tail, "_", score_percentile, ".tsv"))
      write_tsv(signi_drug_df, out_file)
    }
  }

  return(list(
    score_df_list = signi_drugs_score_df_all,
    score_drug_list = signi_drugs_score_list_all
  ))
}

summarize_aggrSig_from_extraction_result <- function(res,
                                                     base_dir = "results",
                                                     metadata_path = here("data/metadata/repurposing_drugs_20200324.csv")) {
  require(dplyr)
  require(readr)
  require(here)

  if (!("score_drug_list" %in% names(res))) {
    stop("Input object must contain 'score_drug_list' from extract_aggrSig_significant_drugs_by_tech()")
  }

  score_list <- res$score_drug_list
  all_keys <- names(score_list)
  technologies <- unique(sub("_.*", "", all_keys))

  # Load metadata once
  drug_info_df <- read.csv(metadata_path, skip = 9)
  results <- list()

  for (tech in technologies) {
    message(paste("Processing:", tech))

    score_keys <- grep(paste0("^", tech, "_"), all_keys, value = TRUE)
    sublist <- score_list[score_keys]

    sub_all_drugs <- unique(unlist(sublist))
    sub_all_drugs_df <- data.frame(significant_drug = sub_all_drugs, stringsAsFactors = FALSE)

    # Add indicator columns for each score
    for (score in names(sublist)) {
      method <- sub(paste0(tech, "_"), "", score)
      sub_all_drugs_df[[method]] <- as.integer(sub_all_drugs_df$significant_drug %in% sublist[[score]])
    }

    # Summarize method groups
    sub_all_drugs_df$LINCS <- rowSums(sub_all_drugs_df[, c("WCS", "NCS", "Tau")], na.rm = TRUE)
    sub_all_drugs_df$Cor <- rowSums(sub_all_drugs_df[, c("Cor_spearman", "Cor_pearson")], na.rm = TRUE)
    sub_all_drugs_df$occurrence <- rowSums(sub_all_drugs_df[, c("CMAP", "LINCS", "Cor")], na.rm = TRUE)

    summary_df <- sub_all_drugs_df[, c("significant_drug", "CMAP", "LINCS", "Cor", "occurrence")]

    # Merge with metadata
    merged_df <- merge(summary_df, drug_info_df,
                       by.x = "significant_drug", by.y = "pert_iname",
                       all.x = TRUE)
    merged_df <- merged_df[order(-merged_df$occurrence), ]
    rownames(merged_df) <- NULL

    # Filter top drugs supported by ≥2 methods
    top_df <- merged_df[rowSums(merged_df[, c("CMAP", "LINCS", "Cor")] > 0, na.rm = TRUE) >= 2, ]
    rownames(top_df) <- NULL

    # Save to disk
    output_dir <- here(base_dir,tech, "03_methodwise")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    write_tsv(merged_df, file.path(output_dir, paste0(tech,"_aggr_drug_summary.tsv")))
    write_tsv(top_df, file.path(output_dir, paste0(tech,"_aggr_top_drugs.tsv")))

    results[[tech]] <- list(all_summary = merged_df, top_drugs = top_df)
  }

  return(results)
}