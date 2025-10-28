# script to perform differential gene expression on processed microarray data matrices
# using limma package: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
# ref (design matrix): https://rpubs.com/ge600/limma
# ref (microarray DE analysis on 2 groups): https://alexslemonade.github.io/refinebio-examples/02-microarray/differential-expression_microarray_01_2-groups.html
# last modified: 07/15/25 - LT
# author: Kewalin Samart

library(limma)
library(magrittr)
library(dplyr)
library(readr)
library(stringr)
library(here)

## set up arguments
args <- commandArgs(TRUE)

# args[1] <- "data/v2/microarray_data_forDE/clean_TB_sample_metadata_classification.tsv"
# args[2] <- 0.05

# --- Check for args ---
if (length(args) < 2) {
  stop("
  Usage:
    Rscript microarray_DE_with_limma.R <metadata_file.tsv> <padj_cutoff>

  Example:
    Rscript microarray_DE_with_limma.R data/v2/microarray_data_forDE/clean_TB_sample_metadata_classification.tsv 0.05
  ")
}

# read in argument file
meta_class_file_path <- args[1]
padj_cutoff <- ifelse(length(args) >= 2, as.numeric(args[2]), 0.05) # default 0.05
lincs_genes <- read.delim(here("data/metadata/LINCSGeneSpaceSub.txt"), sep = "\t")

signature_boolean <- list()
platform_list <- list()

# 1.  Read the batch table
meta_class_df <- read_tsv(meta_class_file_path, show_col_types = FALSE)
study_df <- meta_class_df %>% dplyr::distinct(series_id, SIGNATURE_NAME, EXPRMAT_PATH)
signature_boolean <- logical(nrow(study_df))
up_genes_num <- integer(nrow(study_df))
dn_genes_num <- integer(nrow(study_df))

for (i in 1:nrow(study_df)) {
  # define tag i.e., study id
  tag <- study_df$series_id[i]

  # get the current SIGNATURE_NAME for this loop iteration
  current_signature <- study_df$SIGNATURE_NAME[i]
  message("\n== Processing: ", tag, " | Signature: ", current_signature, " ==")

  # read in expression matrix
  exprmat_path <- paste0(study_df$EXPRMAT_PATH[i])
  expr_mat <- read.delim(here(exprmat_path), sep = "\t")

  # standardise GSM IDs: trim whitespace and set to upper‑case
  colnames(expr_mat) <- trimws(toupper(colnames(expr_mat)))
  meta_class_df$series_id <- trimws(toupper(meta_class_df$series_id))

  # ---- 3.1b  Pull in the classification column from args_df ----
  # Filter for *both* series_id AND the current_signature
  class_df <- meta_class_df %>%
    filter(series_id == tag, SIGNATURE_NAME == current_signature) %>% # filter for samples by signature
    mutate(geo_accession = trimws(toupper(geo_accession))) # normalise case/space

  # here we need a column for condition so that we can focus on case and control in the next step
  class_df <- dplyr::select(class_df, geo_accession, platform_id, CLASSIFICATION) # using the classification column here

  # ---- 3.2  Attach classification and ensure matched GSMs ----
  # Keep only GSMs present in *both* metadata and matrix
  common_ids_infected <- intersect(class_df[trimws(tolower(class_df$CLASSIFICATION)) == "disease without treatment", ]$geo_accession, colnames(expr_mat))
  common_ids_control <- intersect(class_df[trimws(tolower(class_df$CLASSIFICATION)) == "healthy control without treatment", ]$geo_accession, colnames(expr_mat))

  if ((length(common_ids_infected) < 3) | (length(common_ids_control) < 3)) {
    warning("  skipped: either infected or control condition has < 3 classified samples present in the count matrix")
    signature_boolean[i] <- FALSE
    next
  }

  # get all the common ids
  common_ids <- c(common_ids_infected, common_ids_control)
  class_df <- class_df[match(common_ids, class_df$geo_accession), ]
  expr_mat <- expr_mat[, common_ids, drop = FALSE]

  # ---- 3.3  Build the condition factor (now classification exists!) ----
  # get infected/control sample indices
  inf_indices <- which(trimws(tolower(class_df$CLASSIFICATION)) == "disease without treatment")
  control_indices <- which(trimws(tolower(class_df$CLASSIFICATION)) == "healthy control without treatment")

  # Create classification labels
  conditions <- c()
  conditions[inf_indices] <- "infected"
  conditions[control_indices] <- "control"

  # Create a factor for group assignment
  class_df$CLASSIFICATION <- factor(conditions, levels = c("infected", "control"))

  # Platform batch info
  platform_ids <- unique(class_df$platform_id)

  if (length(platform_ids) > 1) {
    print("removing batch effect caused by platform")
    message("Platform IDs: ", paste(platform_ids, collapse = "_"))

    # Ensure platform_id is a factor
    class_df$platform_id <- factor(class_df$platform_id)

    # Use CLASSIFICATION and platform_id to build design matrix
    design_mat <- model.matrix(~ 0 + CLASSIFICATION + platform_id, data = class_df)
  } else {
    design_mat <- model.matrix(~ 0 + CLASSIFICATION, data = class_df)
    print("no batch effect caused by GSE or platform")
  }

  # Rename classification columns only
  group_levels <- levels(class_df$CLASSIFICATION)
  group_cols <- grep("^CLASSIFICATION", colnames(design_mat))
  colnames(design_mat)[group_cols] <- tolower(group_levels)

  ## perform differential gene expression analysis
  # apply linear model to the expression matrix
  fit <- lmFit(expr_mat, design = design_mat)

  # Create the contrast matrix to compare infected vs control
  cont_mat <- makeContrasts(infected_vs_control = infected - control, levels = design_mat)
  fit2 <- contrasts.fit(fit, cont_mat)
  # Apply empirical Bayes to smooth standard errors on the contrasted fit
  fit2 <- eBayes(fit2)
  # Obtain result table for the specified contrast.
  stats_df <- topTable(fit2, coef = "infected_vs_control", number = nrow(expr_mat)) %>%
    tibble::rownames_to_column("Gene")

  ## add gene annotations
  # get gene annotation table
  annot_path <- paste0(here("data/metadata/Homo_sapiens.gene_info.csv"))
  annotdata <- read.delim(annot_path, sep = ",") # GeneID, Symbol, Ensembl
  entrezids <- stats_df$Gene
  annotdata_subset <- annotdata %>% filter(as.character(annotdata$GeneID) %in% entrezids)
  annotdata_subset <- annotdata_subset[, c("GeneID", "Symbol", "Ensembl")]
  res_df <- merge(stats_df, annotdata_subset, by.x = "Gene", by.y = "GeneID", all.x = TRUE, all.y = FALSE)

  # add mean expression values
  group_means <- as.data.frame(fit$coefficients) %>%
    dplyr::select(any_of(c("infected", "control"))) %>%
    tibble::rownames_to_column("Gene")
  res_df <- merge(res_df, group_means, by = "Gene")

  # Select and rename columns for the final results table.
  res_df <- res_df[, c("Gene", "Symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "infected", "control")]
  colnames(res_df)[1] <- "GeneID"
  colnames(res_df)[3] <- "log2FoldChange"

  # identify duplicated genes (EntrezID)
  duplicates_booleans <- duplicated(res_df$GeneID)
  duplicated_genes <- res_df$GeneID[duplicates_booleans]
  duplicated_genes <- duplicated_genes[!is.na(duplicated_genes)]

  if (length(duplicated_genes) != 0) {
    for (gene in duplicated_genes) {
      # look for its maximum logFC value and keep the max
      gene_indices <- which(res_df$GeneID == gene)
      logFC_vec <- res_df$log2FoldChange[gene_indices]
      names(logFC_vec) <- gene_indices
      max_index <- as.numeric(names(logFC_vec)[which.max(abs(logFC_vec))])
      not_max_index <- gene_indices[which(gene_indices != max_index)]
      res_df <- res_df[-c(not_max_index), ]
    }
  }

  # remove NA logFC
  res_df <- res_df[!is.na(res_df$log2FoldChange), ]

  # ---- 3.6  Landmark genes filter ----
  landmark_genes <- as.character(lincs_genes[lincs_genes$Type == "landmark", ]$Entrez.ID)
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
  } else if (nrow(up_df) == 0) {
    message("  No significant up genes (padj < ", padj_cutoff, ")")
    signature_boolean[i] <- FALSE
  }

  if (nrow(dn_df) > 0) {
    dn_genes_num[i] <- nrow(dn_df)
  } else if (nrow(dn_df) == 0) {
    message("  No significant dn genes (padj < ", padj_cutoff, ")")
    signature_boolean[i] <- FALSE
  }

  # save the signatures
  today <- format(Sys.Date(), "%Y%m%d")
  base_fname <- study_df$SIGNATURE_NAME[i]

  # Create output directories (added this for safety)
  dir.create(here("data/v2/DE_results/microarray"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("data/v2/signatures/microarray/full"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("data/v2/signatures/microarray/up"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("data/v2/signatures/microarray/dn"), recursive = TRUE, showWarnings = FALSE)


  if (nrow(res_df) > 0) {
    readr::write_tsv(
      res_df %>%
        dplyr::arrange(dplyr::desc(log2FoldChange)),
      file.path(here::here("data/v2/DE_results/microarray"), paste0(base_fname, "_limma.tsv"))
    )
  }
  if (nrow(sig_df) > 0) {
    readr::write_tsv(
      sig_df %>% arrange(desc(log2FoldChange)),
      file.path(here("data/v2/signatures/microarray/full"), paste0(base_fname, "_full.tsv"))
    )
  }
  if (nrow(up_df) > 0) {
    readr::write_tsv(
      up_df %>% arrange(desc(log2FoldChange)),
      file.path(here("data/v2/signatures/microarray/up"), paste0(base_fname, "_up.tsv"))
    )
  }
  if (nrow(dn_df) > 0) {
    readr::write_tsv(
      dn_df %>% arrange(desc(log2FoldChange)),
      file.path(here("data/v2/signatures/microarray/dn"), paste0(base_fname, "_dn.tsv"))
    )
  }

  message(
    "  ✔  Saved ", nrow(up_df), " up‑regulated and ",
    nrow(dn_df), " down‑regulated genes"
  )

  signature_boolean[i] <- TRUE
}

# 4.  Run summary
study_df$signature <- as.integer(signature_boolean)
study_df$up_genes_num <- as.integer(up_genes_num)
study_df$dn_genes_num <- as.integer(dn_genes_num)

run_info <- file.path(
  here("data/v2/signatures"),
  paste0("microarray_TB_signature_run_info.tsv")
)
readr::write_tsv(study_df, run_info)
message("\nFinished.  Summary written to ", run_info)
