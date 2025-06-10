# script to compute a weighted average gene vector -- an aggregated signature
# created date: 08/31/21
# last modified: 10/04/23
# Kewalin Samart

# load needed
library(readr)
source("./scripts/01_signature_aggregation_functions.R")

# set up arguments
args <- commandArgs(TRUE)
metadata_path <- args[1] # "/data/scratch/samartk/drugrep_tb/inputs/TB_microarray_args.tsv"
data_path <- args[2] # "/data/scratch/samartk/data/uniformly_processed/microarray/signatures/"
direction <- args[3] # direction of differential genes: "up", "dn", or "full"
bg_source <- args[4] # name of the source for background genes to use: "LINCS", "KEGG", "GO", "input data"
output_dir <- args[5] # path to the output directory
extra_arg <- args[6] # extra argument; optional if "LINCS" is specified for bg_source; could be one of the options below or a combination of them as a single string with comma:
                     # (i) "landmark" (by default) (ii) "inferred" (iii) "best inferred" (iv) "not inferred" (v) "reference" for examples: "landmark" or "landmark,inffered,best inferred"
                     # if "input data" is specified for bg_source, this arg could be one of the followings: "up", "dn", "full", and "" (by default)
threshold <- args[7] # aggregated gene score cutoff; 0.4 by default
if(threshold == None){
    threshold <- 0.4
}

membership_matrix = compute_membership_matrix(metadata_path, data_path, direction, bg_source, output_dir, extra_arg)
jaccard_matrix = compute_jaccard_matrix(metadata_path, data_path,direction, output_dir)
aggr_signature = aggregate_signatures(membership_matrix, jaccard_matrix, output_dir, threshold=threshold)
