# multids_RNAseq_DEwithDESeq2.R
# ------------------------------------------------------------
# Batch differential‑expression analysis for multiple RNA‑seq
# datasets using DESeq2.
#
# The per‑dataset *metadata* already contains a cleaned
# classification column with exactly two possible values:
#   1) healthy control without treatment
#   2) disease without treatment
#
# The contrast within a study:
#     disease_without_treatment  vs  healthy_control_without_treatment
# with unique sample conditions
#
# Required R packages: DESeq2, readr, dplyr, AnnotationDbi, org.Hs.eg.db
#
# --------------------------------------------------------------------
# Usage:
#   Rscript multids_RNAseq_DEwithDESeq2.R <args_file.tsv> [padj_cutoff]
#
#   meta_class_file.tsv  : Tab‑separated file listing all datasets.
#                    Mandatory columns:
#                       series_id       (GEO study identifier)
#                       geo_accession   (GEO sample identifier )
#                       SIGNATURE_NAME  (name of signature containing unique sample conditions)
#                       EXPRMAT_PATH    (path to raw‑count matrix TSV)
#                       CLASSIFICATION  (labels:  disease_without_treatment  or  healthy_control_without_treatment)
#   padj_cutoff    : Adjusted‑p significance threshold (default 0.05)
# --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(here)
})

# --------------------------------------------------------------------
# 0.  Parse command‑line arguments
# --------------------------------------------------------------------
argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) < 1) {
  stop("
  Usage:
    Rscript multids_RNAseq_DEwithDESeq2.R <meta_class_file.tsv> [padj_cutoff]

  See header comments for required columns in meta_class_file.tsv.
  ")
}

meta_class_file_path <- argv[1]
padj_cutoff <- ifelse(length(argv) >= 2, as.numeric(argv[2]), 0.05)

######### ------------ Example arguments ------------#########

# meta_class_file_path <- "./data/RNAseq_data_forDE/clean_TB_sample_metadata_classification.tsv"
# padj_cutoff <- 0.05

# 1.  Read the batch table
meta_class_df <- read_tsv(meta_class_file_path, show_col_types = FALSE)
study_df <- meta_class_df %>% dplyr::distinct(series_id, SIGNATURE_NAME, EXPRMAT_PATH)
signature_boolean <- logical(nrow(study_df))
up_genes_num <- integer(nrow(study_df))
dn_genes_num <- integer(nrow(study_df))

# 2.  Helper to coerce expression table to integer matrix
make_expression_matrix <- function(tbl) {
  rn <- tbl[[1]]
  mat <- as.matrix(tbl[-1])
  rownames(mat) <- rn
  storage.mode(mat) <- "integer"
  mat
}

# global variables
landmark_genes_df <- read_tsv(here("data/metadata/LINCSGeneSpaceSub.txt"))
landmark_genes <- as.character(landmark_genes_df[landmark_genes_df$Type == "landmark", ]$`Entrez ID`)


# 3.  Main loop
for (i in seq_len(nrow(study_df))) {
  message("\n== Dataset ", i, " / ", nrow(study_df), " ==")

  tag <- study_df$series_id[i]
  expr_path <- study_df$EXPRMAT_PATH[i]

  # get the current signature for this loop iteration
  current_signature <- study_df$SIGNATURE_NAME[i]
  message("  Study: ", tag, " | Signature: ", current_signature)


  # ---- 3.1  Load inputs ----
  # filter metadata for *both* study AND signature
  metadata <- meta_class_df %>%
    filter(series_id == tag, SIGNATURE_NAME == current_signature) %>% # filters samples by signature
    mutate(geo_accession = trimws(toupper(geo_accession))) # Standardize IDs

  expr_tbl <- read_tsv(here(expr_path), show_col_types = FALSE)
  expr_mat <- make_expression_matrix(expr_tbl)

  # Standardise GSM IDs: trim whitespace and set to upper‑case
  colnames(expr_mat) <- trimws(toupper(colnames(expr_mat)))
  # 3.2  Attach classification and ensure matched GSMs
  # Keep only GSMs present in *both* metadata and matrix
  # This line is now correct because 'metadata' is the filtered subset
  common_ids <- intersect(metadata$geo_accession, colnames(expr_mat))

  if (length(common_ids) < 4) {
    warning("  skipped: < 4 classified samples present in the count matrix")
    signature_boolean[i] <- FALSE
    next
  }

  # Trim & reorder metadata and matrix to identical GSM sets
  metadata <- metadata[match(common_ids, metadata$geo_accession), ]
  expr_mat <- expr_mat[, common_ids, drop = FALSE]

  # 3.3  Build the condition factor (now classification exists!)
  metadata <- metadata %>%
    mutate(condition = dplyr::case_when(
      CLASSIFICATION == "disease without treatment" ~ "disease",
      CLASSIFICATION == "healthy control without treatment" ~ "control",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(condition)) %>%
    mutate(condition = factor(condition, levels = c("control", "disease")))

  if (dplyr::n_distinct(metadata$condition) != 2) {
    warning("  skipped: both groups not present after filtering")
    signature_boolean[i] <- FALSE
    next
  }

  # Check for missing rownames
  sum(is.na(rownames(expr_mat))) # Should be 0
  sum(rownames(expr_mat) == "") # Should also be 0

  # Drop rows with NA in rownames
  expr_mat <- expr_mat[!is.na(rownames(expr_mat)), ]

  # Drop duplicated gene names, keeping the first occurrence
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]


  # ---- 3.4  Differential expression ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = expr_mat,
    colData   = metadata,
    design    = ~condition
  )
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("condition", "disease", "control"))

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column(var = "ENSEMBL") %>%
    arrange(padj)

  # ---- 3.5  Annotate with Entrez & Ensembl IDs ----
  res_df <- res_df %>%
    dplyr::mutate(
      EntrezID = AnnotationDbi::mapIds(org.Hs.eg.db,
        keys = ENSEMBL,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
      ),
      EnsemblID = AnnotationDbi::mapIds(org.Hs.eg.db,
        keys = ENSEMBL,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )

  # rename symbol column
  res_df <- res_df %>%
    dplyr::rename(Symbol = EnsemblID) %>%
    dplyr::filter(!is.na(EntrezID))
  # rename to uniform column names
  colnames(res_df) <- c("Ensembl", "baseMean", "log2FoldChange", "lfcSE", "stat", "P.Value", "adj.P.Val", "GeneID", "Symbol")

  # ---- 3.6  Landmark genes filter ----
  sig_df <- dplyr::filter(res_df, GeneID %in% landmark_genes)

  # ---- 3.7  Significance filter ----
  sig_df <- dplyr::filter(sig_df, !is.na(adj.P.Val) & adj.P.Val < padj_cutoff)

  if (nrow(sig_df) == 0) {
    message("  No significant genes (padj < ", padj_cutoff, ")")
    signature_boolean[i] <- FALSE
    next
  }

  # ---- 3.8  Split Up / Down ----
  up_df <- dplyr::filter(sig_df, log2FoldChange > 0)
  dn_df <- dplyr::filter(sig_df, log2FoldChange < 0)

  if (nrow(up_df) > 0) {
    up_genes_num[i] <- nrow(up_df)
  }

  if (nrow(dn_df) > 0) {
    dn_genes_num[i] <- nrow(dn_df)
  }

  # ---- 3.9  Write outputs ----
  dir.create(here("data/DE_results/RNAseq"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("data/signatures/RNAseq/up"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("data/signatures/RNAseq/dn"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("data/signatures/RNAseq/full"), recursive = TRUE, showWarnings = FALSE)

  today <- format(Sys.Date(), "%Y%m%d")
  base_fname <- study_df$SIGNATURE_NAME[i]
  if (nrow(res_df) > 0) {
    readr::write_tsv(
      res_df %>%
        dplyr::arrange(dplyr::desc(log2FoldChange)),
      file.path(here::here("data/DE_results/RNAseq"), paste0(base_fname, "_DESeq2.tsv"))
    )
  }
  if (nrow(sig_df) > 0) {
    readr::write_tsv(
      sig_df %>% arrange(desc(log2FoldChange)),
      file.path(here("data/signatures/RNAseq/full"), paste0(base_fname, "_full.tsv"))
    )
  }
  if (nrow(up_df) > 0) {
    readr::write_tsv(
      up_df %>% arrange(desc(log2FoldChange)),
      file.path(here("data/signatures/RNAseq/up"), paste0(base_fname, "_up.tsv"))
    )
  }
  if (nrow(dn_df) > 0) {
    readr::write_tsv(
      dn_df %>% arrange(desc(log2FoldChange)),
      file.path(here("data/signatures/RNAseq/dn"), paste0(base_fname, "_dn.tsv"))
    )
  }

  message(
    "  ✔  Saved ", nrow(up_df), " up‑regulated and ",
    nrow(dn_df), " down‑regulated genes"
  )

  if (nrow(up_df) > 0 && (nrow(dn_df) > 0)){
    signature_boolean[i] <- TRUE
  }

}

# 4.  Run summary
study_df$signature <- as.integer(signature_boolean)
study_df$up_genes_num <- as.integer(up_genes_num)
study_df$dn_genes_num <- as.integer(dn_genes_num)

run_info <- file.path(
  here("data/signatures"),
  paste0("signature_run_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv")
)
readr::write_tsv(study_df, run_info)
message("\nFinished.  Summary written to ", run_info)
