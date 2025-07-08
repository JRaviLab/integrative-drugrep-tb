# functions to finalize drug results from individual disease signatures
# last modified: 07/07/25
# Kewalin Samart

# to be documented
get_drug_results <- function(data_to_run, drug_res_path, score_method, score) {
  #' @description this function obtains all the drug results prioritized by a connectivity scores for a set of individual signatures
  #' @param data_to_run a data detail table containing SIGNATURE_NAME column corresponding to signature file names to read in
  #' @param drug_res_path a path to the drug results (start with a first-level folder name inside the project repository)
  #' @param score_method a string indicating a method category: "CMAP", "LINCS", "Cor"
  #' @param score a string indicating a connectivity score name: "CMAP", "WCS", "NCS", "Tau", "Cor_pearson", "Cor_spearman"
  #' @returns combined_drug_df: a dataframe combining all srug results across individual signatures and the six connectivity scores
  #' @author Kewalin Samart
  require(stringr)

  for (i in 1:nrow(data_to_run)) {
    filename <- data_to_run$SIGNATURE_NAME[i]

    # Use score_method as file prefix for LINCS scores (WCS, NCS, Tau), else use score
    file_prefix <- if (score_method == "LINCS") "LINCS" else score

    file_path <- file.path(drug_res_path, paste0(file_prefix, "_", filename, ".tsv"))

    # Read the drug result file
    if (!file.exists(here::here(file_path))) {
      warning(paste("File not found:", file_path))
      next
    }

    drug_df <- read.delim(here::here(file_path), sep = "\t")

    # Extract correct score column based on the actual score value
    if (score == "CMAP") {
      drug_df <- drug_df[, c("scaled_score", "pert", "cell", "t_gn_sym")]
    } else if (score == "WCS") {
      drug_df <- drug_df[, c("WTCS", "pert", "cell", "t_gn_sym")]
    } else if (score == "NCS") {
      drug_df <- drug_df[, c("NCS", "pert", "cell", "t_gn_sym")]
    } else if (score == "Tau") {
      drug_df <- drug_df[, c("Tau", "pert", "cell", "t_gn_sym")]
    } else if (score %in% c("Cor_spearman", "Cor_pearson")) {
      drug_df <- drug_df[, c("cor_score", "pert", "cell", "t_gn_sym")]
    }

    # Add pert_cell and GSE_platform columns
    drug_df$pert_cell <- str_c(drug_df$pert, "_", drug_df$cell)
    drug_df$GSE_platform <- filename

    if (i == 1) {
      combined_drug_df <- drug_df
    } else {
      combined_drug_df <- rbind(combined_drug_df, drug_df)
    }
  }

  return(combined_drug_df)
}


get_DrugDis_ScoreMatrix <- function(combined_drug_df, score, stats = "min"){
  #' @description this function converts the combined_drug_df to a dataframe (matrix filled with scores) where rows: drugs and cols: disease signatures
  #' @param combined_drug_df a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @param score a string indicating a string indicating a choice of connectivity scores: "CMAP","WCS","NCS","Tau","Cor_spearman", "Cor_pearson"
  #' @param stats a string indicating the stats of score to include in the matrix: "median", "min", max".; "min" by default.
  #' @returns transformed_combined_drug_df: a dataframe with rows: drugs, columns: disease signatures, entries: scores
  #' @author Kewalin Samart

  # modify the dataframe for filtering significant drugs
  unique_pert <- unique(combined_drug_df$pert)
  unique_signature <- unique(combined_drug_df$GSE_platform)

  transformed_combined_drug_df <- data.frame(unique_pert = unique_pert, stringsAsFactors = FALSE)

  col_idx = 1
  for(sig in unique_signature){
    print(sig)
    col_idx = col_idx + 1
    # get tau vector for each signature
    score_vec = c()
    for(pert in unique_pert){

      pert_df = combined_drug_df[which(combined_drug_df$pert == pert & combined_drug_df$GSE_platform == sig),]

      if (nrow(pert_df) == 0) {
        score_vec <- append(score_vec, NA)
        next
      }

      if(score == "Tau"){
        if(stats == "min"){
          index = which.min(pert_df$Tau)
          score_vec <- append(score_vec, pert_df$Tau[index])
        }else if(stats == "median"){
          if (nrow(pert_df) == 0 || all(is.na(pert_df$Tau))) {
          score_vec <- append(score_vec, NA)
        } else {
          score_vec <- append(score_vec, median(pert_df$Tau, na.rm = TRUE))
        }

        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$Tau))
        }
      }
      if(score == "NCS"){
        if(stats == "min"){
          index = which.min(pert_df$NCS)
          score_vec <- append(score_vec, pert_df$NCS[index])
        }else if(stats == "median"){
          if (nrow(pert_df) == 0 || all(is.na(pert_df$NCS))) {
            score_vec <- append(score_vec, NA)
          } else {
            score_vec <- append(score_vec, median(pert_df$NCS, na.rm = TRUE))
          }
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$NCS))
        }
      }
      if(score == "WCS"){
        if(stats == "min"){
          index = which.min(pert_df$WTCS)
          score_vec <- append(score_vec, pert_df$WTCS[index])
        }else if(stats == "median"){
          if (nrow(pert_df) == 0 || all(is.na(pert_df$WTCS))) {
            score_vec <- append(score_vec, NA)
          } else {
            score_vec <- append(score_vec, median(pert_df$WTCS, na.rm = TRUE))
          }
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$WTCS))
        }
      }
      if(score == "CMAP"){
        if(stats == "min"){
          index = which.min(pert_df$scaled_score)
          score_vec <- append(score_vec, pert_df$scaled_score[index])
        }else if(stats == "median"){
          if (nrow(pert_df) == 0 || all(is.na(pert_df$scaled_score))) {
            score_vec <- append(score_vec, NA)
          } else {
            score_vec <- append(score_vec, median(pert_df$scaled_score, na.rm = TRUE))
          }
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$scaled_score))
        }
      }
      if(score == "Cor_spearman" | score == "Cor_pearson"){
        if(stats == "min"){
          index = which.min(pert_df$cor_score)
          score_vec <- append(score_vec, pert_df$cor_score[index])
        }else if(stats == "median"){
          if (nrow(pert_df) == 0 || all(is.na(pert_df$cor_score))) {
            score_vec <- append(score_vec, NA)
          } else {
            score_vec <- append(score_vec, median(pert_df$cor_score, na.rm = TRUE))
          }
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$cor_score))
        }
      }
    }
    # now we have the vector; we want to add this as a new column to the data frame
    colnames(transformed_combined_drug_df)[1] <- "unique_pert"
    transformed_combined_drug_df[col_idx] <- score_vec
    colnames(transformed_combined_drug_df)[col_idx] <- sig
  }
  print(transformed_combined_drug_df)

  return(transformed_combined_drug_df)
}

get_DrugCellDis_ScoreMatrix <- function(combined_drug_df, score){
  #' @description this function converts the combined_drug_df to a dataframe (matrix filled with scores) where rows: drug-cell combinations and cols: disease signatures
  #' @param combined_drug_df a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @param score a string indicating a string indicating a choice of connectivity scores: "CMAP","WCS","NCS","Tau","Cor_spearman", "Cor_pearson"
  #' @returns transformed_combined_drug_df: a dataframe with rows: drug-cell combinations, columns: disease signatures, entries: scores
  #' @author Kewalin Samart

  # modify the dataframe for filtering significant drugs; considering cell types
  unique_pert_cell <- unique(combined_drug_df$pert_cell)
  unique_signature <- unique(combined_drug_df$GSE_platform)

  transformed_combined_drugcell_mat = matrix(NA, nrow = length(unique_pert_cell), ncol = length(unique_signature),
                                             dimnames = list(unique_pert_cell, unique_signature))
  col_idx = 1
  for(sig in unique_signature){
    row_idx = 1
    for(pert_cell in unique_pert_cell){
      pert_cell_df = combined_drug_df[which(combined_drug_df$pert_cell == pert_cell & combined_drug_df$GSE_platform == sig),]
      if(nrow(pert_cell_df) != 0){
        if(score == "Tau"){
          score_val = pert_cell_df$Tau[1]
        }
        if(score == "NCS"){
          score_val = pert_cell_df$NCS[1]
        }
        if(score == "WCS"){
          score_val = pert_cell_df$WTCS[1]
        }
        if(score == "CMAP"){
          score_val = pert_cell_df$scaled_score[1]
        }
        if(score == "Cor_spearman" | score == "Cor_pearson"){
          score_val = pert_cell_df$cor_score[1]
        }
        transformed_combined_drugcell_mat[row_idx,col_idx] = score_val
      }else{
        transformed_combined_drugcell_mat[row_idx,col_idx] = 0
      }
      row_idx = row_idx + 1
    }
    col_idx = col_idx + 1
  }
  transformed_combined_drugcell_df <- as.data.frame(transformed_combined_drugcell_mat)
  transformed_combined_drugcell_df$unique_pert <- row.names(transformed_combined_drugcell_df)
  transformed_combined_drugcell_df <- transformed_combined_drugcell_df[c("unique_pert",unique_signature)]
  row.names(transformed_combined_drugcell_df) <- NULL

  return(transformed_combined_drugcell_df)
}

determine_threshold <- function(score_vec, values, tail,  score_percentile=0.90){
  #' @description this function computes a threshold for getting top drugs by a predefined score percentile
  #' @param score_vec a numeric vector of disease-drug scores
  #' @param score_percentile a numeric indicating percentile of choice e.g. 0.75 representing a threshold being the 75th percentile of the flipped score distribution (dist x-axis: positive <---> negative, as negative score preferred for disease-drug reversal)
  #' @param values a string indicating which non-zero values to include in the distribution: "neg", "pos", if other values specified for this argument, it means including all values.
  #' @param tail a string indicating which side of the distribution to get the indicated top percentile from
  #' @returns threshold
  #' @author Kewalin Samart

  # set scores of zero to NA as they are meaningless
  score_vec[score_vec == 0] <- NA
  # considerations of what non-zero values to include in threshold calculation
  if(values == "neg"){
    print("Only consider negative values")
    score_vec[score_vec > 0] <- NA # set positive values to NAs so they would be removed later
  }else if(values == "pos"){
    print("Only consider positive values")
    score_vec[score_vec < 0] <- NA # set positive values to NAs so they would be removed later
  }

  score_vec <- score_vec[!is.na(score_vec)] # remove NAs

  # compute threshold based on a given percentile and tail direction
  if(tail == "left"){
    threshold <- round(unname(quantile(score_vec, 1-score_percentile)[paste0((1-score_percentile)*100,"%")])[1],3)
  }else if(tail == "right"){
    threshold <- round(unname(quantile(score_vec, score_percentile)[paste0((score_percentile)*100,"%")])[1],3)
  }

  print(paste0("threshold:",threshold))

  return(threshold)
}

get_reversing_drugs_freq <- function(transformed_combined_drug_df, threshold, values){
  #' @description this function counts how often the drug got prioritized by the input disease signatures based on a given score threshold
  #' @param transformed_combined_drug_df a dataframe with rows: drugs, columns: disease signatures, entries: scores
  #' @param threshold a negative numeric indicating the threshold for reversal
  #' @param values a string indicating which non-zero values to include in the distribution: "neg", "pos", if other values specified for this argument, it means including all values.
  #' @returns pert_occurrence_df: a sorted dataframe with drug names and their occurrences across the input signatures
  #' @author Kewalin Samart

  # get unique drug names and initialize a dataframe for output
  unique_pert <- unique(transformed_combined_drug_df$unique_pert)
  pert_occurrence_df <- as.data.frame(unique_pert)

  # prepare a matrix-like dataframe for counting occurrences
  rownames(transformed_combined_drug_df) <- transformed_combined_drug_df$unique_pert
  transformed_combined_drug_df$unique_pert <- NULL

  # loop through each drug and disease signatures to get occurrence counts
  pert_above_occurrence <- list()
  for(row in 1:nrow(transformed_combined_drug_df)){
    above_count <- 0
    for(col in 1:ncol(transformed_combined_drug_df)){
      if(values == 'neg'){
        if(transformed_combined_drug_df[row,col] < threshold){
          above_count = above_count + 1
        }
      }else if(values == 'pos'){
        if(transformed_combined_drug_df[row,col] > threshold){
          above_count = above_count + 1
        }
      }else if(values == 'zero'){
        if(abs(transformed_combined_drug_df[row,col]) < threshold){
          above_count = above_count + 1
        }
      }
    }
    pert_above_occurrence <- append(pert_above_occurrence, above_count)
  }
  # add occurrence column to the dataframe
  pert_occurrence_df$occurrence <- as.numeric(pert_above_occurrence)

  # add drug info to the occurrence dataframe
  # "moa","disease_area","t_gn_sym","target","indication"
  # read in drug info from drug repurposing hub (treat this as a extdata)
  drug_metadata <- read.csv(file=here("data/metadata/repurposing_drugs_20200324.csv"),skip = 9)
  pert_occurrence_df <- merge(x=pert_occurrence_df, y=drug_metadata, by.x = "unique_pert", by.y = "pert_iname", all.x=TRUE, all.y=FALSE)
  # select drugs with at least 1 occurrence
  pert_occurrence_df <- pert_occurrence_df[pert_occurrence_df$occurrence > 0,]
  # sort drugs by their occurrence in decreasing order
  pert_occurrence_df = pert_occurrence_df[order(pert_occurrence_df$occurrence, decreasing = TRUE),]

  return(pert_occurrence_df)
}


get_significant_drugs <- function(pert_occurrence_df, transformed_combined_drug_df, percent_reverse = NA, n = NA){
  #' @description Identifies significant drugs based on either a reversal percentage threshold or the top n most frequent occurrences (including ties).
  #' @param pert_occurrence_df A sorted dataframe with drug names and their occurrences.
  #' @param transformed_combined_drug_df A dataframe with drug vs. disease signature scores.
  #' @param percent_reverse A float in [0,1] indicating the percentage of signatures reversed.
  #' @param n An integer: the number of top occurrences to include (ties allowed).
  #' @returns signi_drug_df: A dataframe with significant drugs.

  # check input
  if (!(is.numeric(percent_reverse) | is.numeric(n))) {
    print("Please enter a numerical value for either percent_reverse or n")
    return(NULL)
  }

  # select based on percent_reverse
  if (is.na(n)) {
    print("Selecting drugs by percent_reverse")
    threshold <- percent_reverse * ncol(transformed_combined_drug_df)
    signi_drug_df <- pert_occurrence_df[pert_occurrence_df$occurrence >= threshold, ]
  }
  # select based on top n (including ties at the nth level)
  if (is.na(percent_reverse)) {
    print("Selecting drugs by top n occurrences (including ties)")

    # sort by occurrence
    pert_occurrence_df <- pert_occurrence_df[order(pert_occurrence_df$occurrence, decreasing = TRUE), ]

    # identify the nth highest occurrence value
    if (n <= nrow(pert_occurrence_df)) {
      u <- unique(pert_occurrence_df$occurrence)
      n <- min(n, length(u))
      # Define cutoff: lowest occurrence value to include
      cutoff_occurrence <- u[n]
      print(paste0("occurrence cutoffs: ",cutoff_occurrence))
      signi_drug_df <- pert_occurrence_df[pert_occurrence_df$occurrence >= cutoff_occurrence, ]
    } else {
      signi_drug_df <- pert_occurrence_df
    }
  }

  # check if result is empty
  if (nrow(signi_drug_df) == 0) {
    print("No significant drugs identified. Try adjusting threshold or n.")
  }
  print(paste0("number of significant drugs: ", nrow(signi_drug_df)))
  return(signi_drug_df)
}

get_signi_info <- function(combined_drug_df, signi_drug_df){
  #' @description this function grabs the LINCS signatures information of identified significant drugs from the original combined drug dataframe
  #' @param combined_drug_df a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @param signi_drug_df a dataframe containing the drugs prioritized "significant" based on either the given percent_reverse or (top) n
  #' @returns  signi_info_df: a dataframe with significant drugs and their LINCS signatures information
  #' @author Kewalin Samart
  signi_info_df <- combined_drug_df[which(combined_drug_df$pert %in% signi_drug_df$unique_pert),]

  return(signi_info_df)
}

summarize_significant_drugs <- function(
    metadata_path,
    drug_res_base = "results",
    score_method = "Cor",
    score = "Cor_pearson",
    stats = "median",
    score_percentile = 0.9,
    percent_reverse = NA,
    n = 3,
    values = "neg",
    tail = "left"
) {
  #' @description Summarizes significant drugs from RNA-seq or microarray for a given score
  #' @param metadata_path path to signature detail table starting at a first-level folder within the project repository
  #' @param drug_res_base drug result directory; default is "results"
  #' @param score_method a string indicating a method categories: "CMAP", "LINCS", "Cor"
  #' @param score a string indicating a connectivity score name: "CMAP", "WCS", "NCS", "Tau", "Cor_pearson", "Cor_spearman"
  #' @param percent_reverse A float in [0,1] indicating the percentage of signatures reversed
  #' @param n An integer: the number of top occurrences to include (ties allowed)
  #' @param score_percentile a numeric indicating percentile of choice e.g. 0.75 representing a threshold being the 75th percentile of the flipped score distribution (dist x-axis: positive <---> negative, as negative score preferred for disease-drug reversal)
  #' @param values a string indicating which non-zero values to include in the distribution: "neg", "pos", if other values specified for this argument, it means including all values.
  #' @param tail a string indicating which side of the distribution to get the indicated top percentile from
  #' @returns A named list containing all intermediate and final results

  require(readr)

  # load metadata and filter by signature ID
  metadata_df <- read_tsv(here::here(metadata_path), show_col_types = FALSE)
  data_to_run <- metadata_df[metadata_df$signature == 1, ]

  # build drug results path
  technology <- ifelse(grepl("RNAseq", metadata_path), "RNAseq", "microarray")
  drug_res_path <- file.path(drug_res_base, technology, score_method)

  # get drug scores
  combined_drug_df <- get_drug_results(data_to_run, drug_res_path, score_method, score)

  score_matrix <- get_DrugDis_ScoreMatrix(combined_drug_df, score, stats = stats)
  score_vector <- as.numeric(unlist(combined_drug_df[1]))

  # compute threshold and occurrences
  threshold <- determine_threshold(score_vector, values = values, tail = tail, score_percentile = score_percentile)

  message(glue::glue("Processing: {technology} | {score_method} | {score} | threshold: {round(threshold, 3)}"))

  pert_occurrence <- get_reversing_drugs_freq(score_matrix, threshold = threshold, values = values)

  # get top significant drugs
  significant_drugs <- get_significant_drugs(pert_occurrence, score_matrix, percent_reverse = percent_reverse, n = n)

  # return all intermediate objects
  return(list(
    combined_drug_df = combined_drug_df,
    score_matrix = score_matrix,
    threshold = threshold,
    pert_occurrence = pert_occurrence,
    significant_drugs = significant_drugs
  ))
}

create_drug_method_table <- function(metadata_path, stats = "min", score_percentile = 0.8, percent_reverse = NA, n = 3){
  #' @description Build a table of significant drugs from RNA-seq or microarray for a given score with occurrences of reversal across scores
  #' @param metadata_path path to signature detail table starting at a first-level folder within the project repository
  #' @param stats a string indicating the stats of score to include in the matrix: "median", "min", max".; "min" by default.
  #' @param score_percentile a numeric indicating percentile of choice e.g. 0.75 representing a threshold being the 75th percentile of the flipped score distribution (dist x-axis: positive <---> negative, as negative score preferred for disease-drug reversal)
  #' @param percent_reverse A float in [0,1] indicating the percentage of signatures reversed
  #' @param n An integer: the number of top occurrences to include (ties allowed)
  #' @returns all_signi_drugs_df: a table of significant drugs from RNA-seq or microarray for a given score with occurrences of reversal across scores

  # define input metadata and parameters
  scores_list <- list(
    CMAP = "CMAP",
    WCS = "WCS",
    NCS = "NCS",
    Tau = "Tau",
    Cor_pearson = "Cor_pearson",
    Cor_spearman = "Cor_spearman"
  )

  score_method_map <- list(
    CMAP = "CMAP",
    WCS = "LINCS",
    NCS = "LINCS",
    Tau = "LINCS",
    Cor_pearson = "Cor",
    Cor_spearman = "Cor"
  )

  # collect all significant drugs into a named list
  all_signi <- list()

  for (score in names(scores_list)) {
    result <- summarize_significant_drugs(
      metadata_path = metadata_path,
      score_method = score_method_map[[score]],
      score = score,
      stats = stats,
      score_percentile = score_percentile,
      percent_reverse = NA,
      n = n
    )
    if (!is.null(result$significant_drugs) && nrow(result$significant_drugs) > 0) {
      all_signi[[score]] <- result$significant_drugs$unique_pert
    }
  }

  # create all_signi_drugs_df: wide-format binary matrix
  all_drugs <- unique(unlist(all_signi))
  all_signi_drugs_df <- data.frame(significant_drug = all_drugs)

  for (score in names(scores_list)) {
    all_signi_drugs_df[[score]] <- as.integer(all_signi_drugs_df$significant_drug %in% all_signi[[score]])
  }
  return(all_signi_drugs_df)
}

summarize_drugs_bymethods <- function(all_signi_drugs_df, technology, dirname = "results") {
  #' @description Summarizes final significant drugs based on CMAP, LINCS (WCS, NCS, Tau), and Cor (Cor_spearman, Cor_pearson)
  #' @param all_signi_drugs_df A dataframe with binary indicators per score and column "significant_drug"
  #' @param technology A string for output file prefix indicating data technology
  #' @param dirname Output directory; default is "results"
  #' @return indiv_drugs_res_top: Dataframe of top drugs in ≥2 of the 3 method groups
  #' @author Kewalin Samart

  require(readr)
  require(here)

  # handle method groups dynamically
  lincs_scores <- intersect(c("WCS", "NCS", "Tau"), colnames(all_signi_drugs_df))
  cor_scores <- intersect(c("Cor_spearman", "Cor_pearson"), colnames(all_signi_drugs_df))
  cmap_score <- "CMAP"

  # add LINCS and Cor group summary columns
  all_signi_drugs_df$LINCS <- if (length(lincs_scores) > 0) {
    rowSums(all_signi_drugs_df[, lincs_scores, drop = FALSE], na.rm = TRUE)
  } else 0

  all_signi_drugs_df$Cor <- if (length(cor_scores) > 0) {
    rowSums(all_signi_drugs_df[, cor_scores, drop = FALSE], na.rm = TRUE)
  } else 0

  if (!("CMAP" %in% colnames(all_signi_drugs_df))) {
    all_signi_drugs_df$CMAP <- 0
  }

  # compute total across 3 methods
  all_signi_drugs_df$occurrence <- rowSums(all_signi_drugs_df[, c("CMAP", "LINCS", "Cor")], na.rm = TRUE)

  # subset summary columns
  summary_df <- all_signi_drugs_df[, c("significant_drug", "CMAP", "LINCS", "Cor", "occurrence")]

  # add metadata (optional)
  drug_info_path <- here::here("data/metadata/repurposing_drugs_20200324.csv")
  if (file.exists(drug_info_path)) {
    drug_info_df <- read.csv(drug_info_path, skip = 9)
    summary_df <- merge(summary_df, drug_info_df,
                        by.x = "significant_drug", by.y = "pert_iname",
                        all.x = TRUE, all.y = FALSE)
  }

  # sort and write full summary
  summary_df <- summary_df[order(-summary_df$occurrence), ]
  rownames(summary_df) <- NULL

  full_summary_path <- here::here(dirname, paste0(technology, "_indiv_drug_summary.tsv"))
  if (!dir.exists(here::here(dirname))) dir.create(here::here(dirname), recursive = TRUE)
  write_tsv(summary_df, file = full_summary_path)

  # select drugs in ≥2 of the 3 methods
  bin_mat <- summary_df[, c("CMAP", "LINCS", "Cor")]
  bin_mat <- as.data.frame(lapply(bin_mat, function(x) as.integer(x > 0)))
  bin_mat$method_freq <- rowSums(bin_mat)

  indiv_drugs_res_top <- summary_df[bin_mat$method_freq >= 2, ]
  rownames(indiv_drugs_res_top) <- NULL

  top_summary_path <- here::here(dirname, paste0(technology, "_indiv_top_drugs.tsv"))
  write_tsv(indiv_drugs_res_top, file = top_summary_path)

  return(indiv_drugs_res_top)
}
