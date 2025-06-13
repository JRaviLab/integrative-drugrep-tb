# Script to compute a weighted average gene vector -- an aggregated signature
# last modified: 06/13/25
# Kewalin Samart

# load needed
library(readr)
library(here)
source(here("scripts/01_signature_aggregation_functions.R"))

# set up arguments
args = commandArgs(TRUE)

metadata_path = args[1] # path to metadata file for the signatures
data_path = args[2] # path to the signatures to be aggregated
direction = args[3] # direction of differential genes: "up", "dn", or "full" (up+dn)
bg_source = args[4] # name of the source for background genes to use: "LINCS" (recommended for disease-drug comparison),
                     # "KEGG", "GO", "input data"
output_dir = args[5] # path to the output directory
extra_arg     <- ifelse(length(args) >= 6, args[6], "") # extra argument for signature file names; "" by default
threshold     <- ifelse(length(args) >= 7, as.numeric(args[7]), 0.4) # aggregated gene score cutoff; 0.4 by default

###### Example arguments ######
# metadata_path = here("data/TB_microarray_args.tsv")
# data_path = here("data/microarray_TBsignatures")
# direction = "up"
# bg_source = "LINCS" # L1000 genes
# extra_arg= "landmark"
# output_dir = here("data/aggregated_signatures")
# threshold = NULL # aggregated gene score cutoff; 0.4 by default

# 1. define membership matrix
membership_matrix = compute_membership_matrix(metadata_path, data_path, direction, bg_source, output_dir, extra_arg, save_result)
# 2. calculate Jaccard similarity matrix
jaccard_matrix = compute_jaccard_matrix(metadata_path, data_path,direction, output_dir, save_result)
# 3. compute aggregated signature
aggr_signature = aggregate_signatures(membership_matrix, jaccard_matrix, output_dir, threshold=threshold, save_result)

