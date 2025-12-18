# functions for aggregating multiple gene signatures
# created date: 09/22/23
# last modified: 11/13/25
# Kewalin Samart

# import bg functions
source(here("scripts/00_background_genes_PA_functions.R"))

jaccard_score <- function(data1, data2) {
  #' @description This functions calculate jaccard similarity score between different two set of data
  #' @reference jaccard function: https://www.statology.org/jaccard-similarity-in-r/
  #' @param data1 set of data no.1; a vector of data1
  #' @param data2 set of data no.2; a vector of data2
  #' @return jaccard jaccard similarity score
  #' @author Kewalin Samart

  intersection <- length(intersect(data1, data2))
  union <- length(data1) + length(data2) - intersection
  jaccard <- intersection / union

  return(jaccard)
}

compute_membership_matrix <- function(metadata_path,
                                      data_path,
                                      direction,
                                      bg_source,
                                      output_dir,
                                      lincs_genetype = NULL,
                                      extension = ".tsv",
                                      tissue_category = NULL,
                                      sample_origin = NULL,
                                      save_result = TRUE) {
  #' @description Given a set of signatures and source name for background genes, this function compyte a membership matrix based on the set of background genes of choice
  #' @param metadata_path metadata of the input signatures
  #' @param data_path path to signature data e.g., "../data/microarray_TBsignatures"
  #' @param direction a string indicating a regulation direction: "up", "dn", "full" (up+dn)
  #' @param bg_source name of the source for background genes to use: "LINCS", "KEGG", "GO", "input data"
  #' @param output_dir path to the output directory
  #' @param lincs_genetype extra argument; optional if "LINCS" is specified for bg_source; could be one of the options below or a combination of them as a single string with comma:
  #' (i) "landmark" (ii) "inferred" (iii) "best inferred" (iv) "not inferred" (v) "reference" for examples: "landmark" or "landmark, inferred, best inferred"
  #' if "input data" is specified for bg_source, this arg could be one of the followings: "up", "dn", "full", and "" (by default)
  #' @param save_result a boolean indicating user's preference whether to automatically save the result
  #' @returns final_sigval resulting gene membership matrix
  #' @author Kewalin Samart

  # import libraries
  require(tidyverse)
  require(dplyr)
  require(readr)

  # input validation
  stopifnot(
    "metadata_path must be a valid file path" = file.exists(metadata_path),
    "data_path must be a valid directory" = dir.exists(data_path),
    "direction must be 'up', 'dn', or 'full'" = tolower(direction) %in% c("up", "dn", "full")
  )

  # get background genes across datasets
  message(paste0("Getting background genes from : ", bg_source))
  bg_genes <- get_bg_genes(
    bg_source = bg_source,
    metadata_path = metadata_path,
    data_path = data_path,
    direction = direction,
    lincs_genetype = lincs_genetype,
    extension = extension
  )
  bg_genes_df <- as.data.frame(bg_genes)
  colnames(bg_genes_df)[1] <- "GeneID"

  print(paste0("Number of background genes: ", length(bg_genes)))

  data_to_run <- read_tsv(metadata_path)
  data_to_run <- data_to_run[data_to_run$signature == 1, ]
  # filter by any specify context
  # tissue category
  if (!is.null(tissue_category)) {
    data_to_run <- data_to_run[data_to_run$tissue_category == tissue_category, ]
  }
  # sample origin
  if (!is.null(sample_origin)) {
    data_to_run <- data_to_run[data_to_run$sample_origin == sample_origin, ]
  }

  for (i in 1:nrow(data_to_run)) {
    print(paste0("iteration ", i))

    # set full variables
    signature_name <- data_to_run$SIGNATURE_NAME[i]
    signature_filename <- paste0(signature_name, "_", direction, ".tsv")
    print(paste0("reading in up and dn genes: ", signature_filename))
    signature_path <- paste0(data_path, "/", direction, "/", signature_filename)

    print(file.exists(signature_path))
    if (file.exists(signature_path)) {
      gene_signature <- read.delim(signature_path, sep = "\t")
    } else {
      next
    }

    clean_gene_signature <- gene_signature[!duplicated(gene_signature$GeneID), ]

    # compute membership vectors with DE vals
    to_merge_logFC <- clean_gene_signature[c("GeneID", "log2FoldChange")]

    sigval_replaced <- merge(bg_genes_df, to_merge_logFC, by = "GeneID", all.x = TRUE)
    sigval_replaced[is.na(sigval_replaced)] <- 0
    colnames(sigval_replaced)[2] <- signature_filename

    if (i == 1) {
      old_sigval <- sigval_replaced
    } else {
      # merge dataframes
      new_sigval <- merge(old_sigval, sigval_replaced, by = "GeneID")
      old_sigval <- new_sigval
      print(paste0("Gene Membership matrix with DE values 's dimension: ", dim(old_sigval)))
    }
  }
  final_sigval <- new_sigval
  col_names <- colnames(final_sigval)
  final_sigval <- final_sigval[, c("GeneID", col_names[!col_names == "GeneID"])]

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # check output dimension
  print(paste0("Number of DE results/signatures expected: ", nrow(data_to_run)))
  print(paste0(
    "Number of columns in membership matrix (excluding GeneID): ",
    ncol(final_sigval[, !colnames(final_sigval) %in% "GeneID"])
  ))
  print(paste0("The expected output dimension matches: ", nrow(data_to_run) == ncol(final_sigval[, !colnames(final_sigval) %in% "GeneID"])))

  # set row names to GeneID
  row.names(final_sigval) <- final_sigval$GeneID
  final_sigval$GeneID <- NULL
  # reorder columns
  final_sigval <- final_sigval[, order(colnames(final_sigval))]
  final_sigval <- as.matrix(final_sigval)

  if (save_result) {
    saveRDS(
      final_sigval,
      file = file.path(output_dir, paste0(direction, "_membership_mat.rds"))
    )

    df_out <- tibble::rownames_to_column(
      as.data.frame(final_sigval),
      var = "GeneID"
    )

    readr::write_tsv(
      df_out,
      file.path(output_dir, paste0(direction, "_membership_mat.tsv"))
    )

    message(
      "The computed Gene Membership matrix (", direction,
      ") was saved at ", output_dir
    )
  }

  return(final_sigval)
}

compute_jaccard_matrix <- function(metadata_path,
                                   data_path,
                                   direction,
                                   output_dir,
                                   tissue_category = NULL,
                                   sample_origin = NULL,
                                   save_result = TRUE) {
  #' @description This function computes a jaccard similarity matrix of a given set of signatures
  #' @param data_to_run metadata of the input signatures
  #' @param data_path e.g., "/data/scratch/samartk/drugrep_tb/data/uniformly_processed/microarray/signatures/up/"
  #' @param direction a string indicating a regulation direction: "up", "dn", "full" (up+dn)
  #' @param output_dir e.g., "/data/scratch/samartk/drugrep_tb/data/uniformly_processed/microarray/signatures/up/"
  #' @returns mat resulting jaccard similarity matrix
  #' @author Kewalin Samart

  # load libraries
  require(readr)

  # input validation
  stopifnot(
    "metadata_path must be a valid file path" = file.exists(metadata_path),
    "data_path must be a valid directory" = dir.exists(data_path),
    "direction must be 'up', 'dn', or 'full'" = tolower(direction) %in% c("up", "dn", "full")
  )

  # get data_to_run df
  data_to_run <- read_tsv(metadata_path)
  data_to_run <- data_to_run[data_to_run$signature == 1, ]
  # filter by any specify context
  # tissue category
  if (!is.null(tissue_category)) {
    data_to_run <- data_to_run[data_to_run$tissue_category == tissue_category, ]
  }
  # sample origin
  if (!is.null(sample_origin)) {
    data_to_run <- data_to_run[data_to_run$sample_origin == sample_origin, ]
  }

  # initialize empty vectors for storage
  vec <- c()
  name_vec <- c()

  n <- 0
  # loop to get gene sets for each dataset
  for (i in 1:(nrow(data_to_run))) {
    # get accession no. and file name for reading in gene signatures
    signature_name <- data_to_run$SIGNATURE_NAME[i]

    signature_filename <- paste0(signature_name, "_", direction, ".tsv")

    signature_path <- paste0(data_path, "/", direction, "/", signature_filename)

    if (file.exists(signature_path)) {
      gene_signature <- read.delim(signature_path, sep = "\t")
      clean_gene_signature <- gene_signature[!duplicated(gene_signature$GeneID), ]
      GeneIDs <- clean_gene_signature$GeneID

      n <- n + 1
      # append vector of genes to the big vector
      vec[[n]] <- GeneIDs
      # create proper name for each gene list
      name <- signature_filename
      # append names to the vectors of names
      name_vec[[n]] <- name
    } else {
      next
    }
  }

  # name the vectors
  names(vec) <- name_vec
  # reorder the vectors of genes by names
  vec <- vec[order(names(vec))]
  # construct matrix
  mat <- matrix(data = NA, nrow = length(vec), ncol = length(vec))

  # assign col/row names
  colnames(mat) <- names(vec)
  rownames(mat) <- names(vec)

  # compute Jaccard score matrix
  for (r in 1:length(vec)) {
    for (c in 1:length(vec)) {
      if (c == r) {
        mat[r, c] <- 1
      } else if (c > r) {
        mat[r, c] <- jaccard_score(vec[[r]], vec[[c]])
      } else {
        mat[r, c] <- jaccard_score(vec[[r]], vec[[c]])
      }
    }
  }
  # add a column of signature names
  # mat <- as.data.frame(mat)
  # mat$names <- names(vec)
  # mat <- mat[, c("names", names(vec))]
  mat[is.na(mat)] <- 0

  n_expected <- nrow(data_to_run)
  n_actual <- length(vec)

  message("Number of signatures expected (metadata): ", n_expected)
  message("Number of signatures loaded: ", n_actual)
  message("Jaccard matrix dimensions: ", nrow(mat), " x ", ncol(mat))

  if (nrow(mat) != ncol(mat)) {
    stop("Jaccard matrix is not square — this should never happen.")
  }

  if (n_actual != n_expected) {
    warning(
      "Not all signatures in metadata were found on disk. ",
      "Expected: ", n_expected,
      ", loaded: ", n_actual
    )
  }

  if (save_result) {
    saveRDS(
      mat,
      file = file.path(output_dir, paste0(direction, "_jaccard_mat.rds"))
    )

    write_tsv(
      tibble::rownames_to_column(as.data.frame(mat), "signature"),
      file = file.path(output_dir, paste0(direction, "_jaccard_mat.tsv"))
    )

    message("Saved Jaccard matrix (", direction, ") to ", output_dir)
  }

  # reorder columns
  mat <- mat[, order(colnames(mat))]


  return(mat)
}

aggregate_signatures <- function(gene_membership_matrix,
                                 jaccard_matrix,
                                 output_dir,
                                 direction,
                                 threshold_pct = 0.90,
                                 save_result = TRUE) {
  #' @description Given a gene membership and a jaccard similarity matrix computed from a set of signatures, this function compute an aggregated signature.
  #' @param gene_membership_matrix gene membership matrix returned by compute_membership_matrix(...)
  #' @param jaccard_matrix jaccard similarity matrix returned by compute_jaccard_matrix(...)
  #' @param output_dir path to the output directory
  #' @param direction a string indicating a regulation direction: "up", "dn", "full" (up+dn)
  #' @param threshold a thershold for selecting a set of significantly aggregated genes; set to 0.9 by default meaning the selected genes are present in at least 40% of the signatures
  #' @returns selected_genes_df final aggregated signature: a list of genes and their aggregated gene scores
  #' @author Kewalin Samart

  # compute average jaccard scores across signatures
  ## specified disease
  jaccard_mean_vec <- colMeans(jaccard_matrix)

  ### (specified disease gene DE score matrix) x (mean jaccard vector/sum(mean jaccard vector))
  aggregated_gene_sig <- as.data.frame(gene_membership_mat %*% (jaccard_mean_vec / sum(jaccard_mean_vec)))
  colnames(aggregated_gene_sig)[1] <- "aggregated_GeneScores"
  aggregated_gene_sig$GeneID <- row.names(gene_membership_df)

  # calculate threshold from the threshold_pct
  threshold <- quantile(abs(aggregated_gene_sig$aggregated_GeneScores), probs = threshold_pct)
  # select genes with logFC > threshold
  selected_genes_df <- aggregated_gene_sig[abs(aggregated_gene_sig$aggregated_GeneScores) > threshold, ]

  if (save_result) {
    saveRDS(selected_genes_df, file = paste0(output_dir, "/", direction, "_aggregated_signature.rds"))
    write_tsv(selected_genes_df, file = paste0(output_dir, "/", direction, "_aggregated_signature.tsv"))
    print(paste0("The final aggregated signature was saved at ", output_dir))
  }

  return(selected_genes_df)
}
