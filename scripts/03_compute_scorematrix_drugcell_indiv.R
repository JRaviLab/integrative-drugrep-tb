# script to compute score matrix of drug-cell combinations cross individual signatures
# created date: 01/12/24
# Kewalin Samart

library(readr)

source("./scripts/03_summarize_drugs_methodswide_functions.R")

# set up arguments for the summarization workflow
args <- commandArgs(TRUE)
drugres_metadata_path <- args[1] # "./inputs/TB_rnaseq_args.tsv"
drug_summ_path <- args[2] # "./inputs/drug_summarization_rnaseq_args.tsv"
prefix <- args[3]  # "uni_rnaseq"
output_dir <- args[4] # "./results/uniformly_processed/RNAseq/"

# read in metadata for disease signatures and drug results
data_to_run <- read.delim(drugres_metadata_path, sep="\t")
data_to_run <- data_to_run[which(data_to_run$signature == 1), ]
print("------------------ Drug Summarization of Individual-Signatures Results ------------------")
print(paste0("Number of drug results to summarize: ",dim(data_to_run)[1]))

drug_summ_args <- read.delim(file=drug_summ_path, sep="\t")

# set up output directory
# make the directory to write the files to (needs to be high I/O capable)
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

for(i in 1:nrow(drug_summ_args)){

  # getting arguments for drug summarization
  drug_res_path <- drug_summ_args$drug_res_path[i]
  score_method <- drug_summ_args$score_method[i]
  score <- drug_summ_args$score[i]
  percent_reverse <- drug_summ_args$percent_reverse[i]
  n <- drug_summ_args$n[i]
  score_percentile <- drug_summ_args$score_percentile[i]

  print(paste0("Getting significant drugs for the disease signatures ", drugres_metadata_path," prioritized by ", score))
  combined_drug_df <- get_drug_results(data_to_run = data_to_run, drug_res_path = drug_res_path, score_method = score_method, score=score)
  transformed_combined_drug_df <- get_DrugCellDis_ScoreMatrix(combined_drug_df,score = score)
  write_tsv(transformed_combined_drug_df,file = paste0(output_dir,"03a_",prefix,"_",score,"_drugcells_indivSig_matrix.tsv"))
}
