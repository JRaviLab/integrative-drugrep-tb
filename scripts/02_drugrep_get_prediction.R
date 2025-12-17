# script to compute connectivity scores and get drug results given input individual disease signatures
# method choices: "LINCS", "CMAP", "Cor_spearman", "Cor_pearson"
# last modified: 11/17/25
# Kewalin Samart
library(here)

# import needed functions
source(here("scripts/01_signature_aggregation_functions.R"))
source(here("scripts/02_signatureSearch_connectivity_scores_functions.R"))

# set up arguments
args <- commandArgs(TRUE)
sig_metadata_path <- args[1] # e.g. "data/signatures/RNASeq_TB_signature_run_info.tsv"
sig_data_path <- args[2] # e.g., "data/signatures/RNAseq"
drugdb_name <- args[3] # "LINCS", "CMAP"
score_method <- args[4] # "LINCS", "CMAP", "Cor_spearman", "Cor_pearson"
output_dir <- args[5] # e.g., "results/RNAseq/LINCS"

sig_metadata_path <- "data/signatures/RNASeq_TB_signature_run_info.tsv"
sig_data_path <- "data/signatures/RNAseq"
drugdb_name <- "LINCS"
score_method <- "LINCS"
output_dir <- "results/RNAseq/LINCS"

# read in metadata file
data_to_run <- read.delim(here(sig_metadata_path), sep = "\t")
data_to_run <- data_to_run[which(data_to_run$signature == 1), ]
print(paste0("Number of signtures to aggregate: ", dim(data_to_run)[1]))

# set up drug database
if (drugdb_name == "CMAP") {
  db_path <- setup_cmap1db()
} else if (drugdb_name %in% c("LINCS")) {
  db_path <- setup_lincsdb()
}

# ## --------- individual signatures ----------- ##
#
# iterate through each signature info
# compare the disease signature against drug profiles using the method of choice
for (i in 1:nrow(data_to_run)) {
  print(paste0("iteration ", i))
  print(paste0("Method: ", score_method))
  # get info from the metadata file
  file_name <- data_to_run$SIGNATURE_NAME[i]
  print(paste0("Quantifying drug candidates for ", file_name))

  if (score_method %in% c("LINCS", "CMAP")) {
    # get signature paths
    up_sig_path <- paste0(sig_data_path, "/up/", file_name, "_up.tsv")
    dn_sig_path <- paste0(sig_data_path, "/dn/", file_name, "_dn.tsv")

    # get disease signatures
    if (file.exists(up_sig_path) && file.exists(dn_sig_path)) { # check if the signatures exist
      # up signature
      up_genes <- get_updn_signature(up_sig_path)
      # dn signature
      dn_genes <- get_updn_signature(dn_sig_path)
    } else {
      next
    }
  } else if (score_method %in% c("Cor_spearman", "Cor_pearson")) {
    full_sig_path <- paste0(sig_data_path, "/full/", file_name, "_full.tsv")
    if (file.exists(full_sig_path)) {
      full_signature_matrix <- get_full_signature(full_sig_path)
    } else {
      next
    }
  }
  # query the input gene signatures against all drug signatures in the selected database
  # finalize the drug result dataframe: add metadata info, filter out non-FDA-approved drugs, and remove rows with NAs
  if (score_method == "CMAP") {
    final_res <- compute_CMap1_scores(up_genes, dn_genes, db_path)
  } else if (score_method == "LINCS") {
    final_res <- compute_CMap2lincs_scores(up_genes, dn_genes, db_path)
  } else if (score_method %in% c("Cor_spearman", "Cor_pearson")) {
    final_res <- compute_Cor_based_scores(full_signature_matrix, score_method, db_path)
  }

  if (!dir.exists(here(output_dir))) {
    dir.create(here(output_dir), recursive = TRUE)
  }

  write_tsv(final_res, file = here(paste0(output_dir, "/", score_method, "_", file_name, ".tsv")))
  print(paste0("Individual signature results saved for ", score_method, "_", file_name))
}
