# script to compute connectivity scores for aggregated signatures
# inputs:
#   1) technology: "microarray" or "RNAseq"
#   2) score_method: "LINCS", "CMAP", "Cor_spearman", "Cor_pearson"
#   3) db_name: "LINCS" or "CMAP"

library(here)
library(readr)

# import functions
source(here("scripts/01_signature_aggregation_functions.R"))
source(here("scripts/02_signatureSearch_connectivity_scores_functions.R"))

# set up arguments (defaults to RNAseq/LINCS if not specified)
args <- commandArgs(TRUE)

tech <- ifelse(length(args) >= 1, args[1], "RNAseq")
score_method <- ifelse(length(args) >= 2, args[2], "LINCS")
db_name <- ifelse(length(args) >= 3, args[3], "LINCS")

print(paste0("Technology: ", tech))
print(paste0("Method: ", score_method))
print(paste0("DB: ", db_name))

# Validate inputs
if (!(tech %in% c("microarray", "RNAseq"))) {
  stop("technology must be: microarray or RNAseq")
}

if (!(score_method %in% c("LINCS", "CMAP", "Cor_spearman", "Cor_pearson"))) {
  stop("score_method must be: LINCS, CMAP, Cor_spearman, Cor_pearson")
}

if (!(db_name %in% c("LINCS", "CMAP"))) {
  stop("db_name must be: LINCS or CMAP")
}

# ## --------- aggregated signatures ----------- ##
# signature and output paths
sig_data_path <- here("data/signatures", tech, "aggregated_signatures")

if (!dir.exists(sig_data_path)) {
  stop("Signature directory does not exist: ", sig_data_path)
}

if (score_method == "CMAP") {
  output_dir <- here("results", tech, "CMAP")
} else if (score_method == "LINCS") {
  output_dir <- here("results", tech, "LINCS")
} else {
  output_dir <- here("results", tech, "Cor")
}

# setup drug database
if (db_name == "CMAP") {
  db_path <- setup_cmap1db()
} else {
  db_path <- setup_lincsdb()
}

# load aggregated signatures
up_sig_path <- file.path(sig_data_path, "up_aggregated_signature.tsv")
dn_sig_path <- file.path(sig_data_path, "dn_aggregated_signature.tsv")

if (!(file.exists(up_sig_path) && file.exists(dn_sig_path))) {
  stop("Aggregated signature files not found at: ", sig_data_path)
}

# prepare inputs
if (score_method %in% c("LINCS", "CMAP")) {

  up_genes <- as.character(read_tsv(up_sig_path, show_col_types = FALSE)$GeneID)
  dn_genes <- as.character(read_tsv(dn_sig_path, show_col_types = FALSE)$GeneID)

} else {

  up_aggr <- read_tsv(up_sig_path, show_col_types = FALSE)
  dn_aggr <- read_tsv(dn_sig_path, show_col_types = FALSE)

  combined <- rbind(up_aggr, dn_aggr)
  combined <- combined[order(-combined$aggregated_GeneScores), ]

  full_signature_matrix <- as.matrix(combined[, !colnames(combined) %in% "GeneID"])
  rownames(full_signature_matrix) <- combined$GeneID
}

# run connectivity scoring
if (score_method == "CMAP") {

  final_res <- compute_CMap1_scores(up_genes, dn_genes, db_path)

} else if (score_method == "LINCS") {

  final_res <- compute_CMap2lincs_scores(up_genes, dn_genes, db_path)

} else {

  final_res <- compute_Cor_based_scores(full_signature_matrix, score_method, db_path)
}

# save prediction output
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

outfile <- file.path(output_dir, paste0(score_method, "_aggregated_signature.tsv"))

write_tsv(final_res, outfile)

print(paste0("Aggregated signature result saved at ", outfile))
