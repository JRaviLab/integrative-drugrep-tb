# script to get drug result summary and top drugs from individual signatures
# created date: 07/20/22
# last modified: 09/08/24
# Kewalin Samart

library(readr)

source("./scripts/03_summarize_drugs_methodswise_functions.R")

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

all_signi_drugs <- c()
signi_drugs_score_list <- list()

for(i in 1:nrow(drug_summ_args)){
  # getting arguments for drug summarization
  drug_res_path <- drug_summ_args$drug_res_path[i]
  score_method <- drug_summ_args$score_method[i]
  score <- drug_summ_args$score[i]
  percent_reverse <- drug_summ_args$percent_reverse[i]
  n <- drug_summ_args$n[i]
  values <- drug_summ_args$values[i]
  tail <- drug_summ_args$tail[i]
  score_percentile <- drug_summ_args$score_percentile[i]

  print(paste0("Getting significant drugs for the disease signatures ", drugres_metadata_path," prioritized by ", score))
  combined_drug_df <- get_drug_results(data_to_run = data_to_run, drug_res_path = drug_res_path, score_method = score_method, score=score)
  score_vec <- as.numeric(unlist(combined_drug_df[1]))
  threshold <- determine_threshold(score_vec, values, tail, score_percentile=score_percentile)

  transformed_combined_drug_df <- get_DrugDis_ScoreMatrix(combined_drug_df,score = score, stats = "max")
  pert_occurrence_df <- get_reversing_drugs_freq(transformed_combined_drug_df, threshold=threshold, values)
  signi_drug_df <- get_significant_drugs(pert_occurrence_df, transformed_combined_drug_df, percent_reverse = percent_reverse, n = n)
  print(paste0("top n:",n))
  #write_tsv(signi_drug_df,file = paste0(output_dir,score,"_top_",score_percentile,"pct_",percent_reverse,"reversed_drugs.tsv"))

  signi_drugs_score <- signi_drug_df$unique_pert
  signi_drugs_score_list[[score]] <- signi_drugs_score # named list containing significant drugs prioritized by different metrics
  all_signi_drugs <- c(all_signi_drugs,signi_drugs_score)

}

all_signi_drugs_df <- data.frame(unique(all_signi_drugs))
colnames(all_signi_drugs_df)[1] <- "significant_drug"

print("Summarizing drug results based on their occurrence of being prioritized by 3 methods: CMAP, LINCS, and Cor methods")

for(score in names(signi_drugs_score_list)){
  score_signi_drugs_boolean <- all_signi_drugs_df$significant_drug %in% signi_drugs_score_list[[score]]
  score_signi_drugs_numeric <- replace(score_signi_drugs_boolean, score_signi_drugs_boolean==TRUE, 1)
  all_signi_drugs_df[score] <- score_signi_drugs_numeric
}

# summarize based on 3 methods: CMAP, LINCS, and Cor methods
all_signi_drugs_df$LINCS = rowSums(all_signi_drugs_df[,c("WCS", "NCS", "Tau")])
all_signi_drugs_df$Cor = rowSums(all_signi_drugs_df[,c("Cor_spearman", "Cor_pearson")])
all_signi_drugs_df$occurrence = rowSums(all_signi_drugs_df[,c("CMAP", "LINCS", "Cor")])
signi_drug_summary_df <- all_signi_drugs_df[,c("significant_drug","CMAP","LINCS","Cor","occurrence")]

drug_info_df <- read.csv(file="./data/metadata/repurposing_drugs_20200324.csv",skip=9)
merged_signi_drug_summary <- merge(signi_drug_summary_df,
                                   drug_info_df,
                                   by.x = "significant_drug",
                                   by.y = "pert_iname",
                                   all.x = TRUE,
                                   all.y = FALSE)
merged_signi_drug_summary <- merged_signi_drug_summary[order(-merged_signi_drug_summary$occurrence),]
rownames(merged_signi_drug_summary) <- 1:nrow(merged_signi_drug_summary)

write_tsv(merged_signi_drug_summary ,file = paste0(output_dir,prefix,"_indiv_drug_summary.tsv"))

# replace non-zero values with 1
indiv_drugs_score <- merged_signi_drug_summary[,c("CMAP","LINCS","Cor")]
indiv_drugs_score[-1] <- as.integer(indiv_drugs_score[-1] != 0)
# subset only drugs that show up in at least 2 out of the 3 methods
indiv_drugs_score$medthod_freq <- rowSums(indiv_drugs_score)
indiv_drugs_res_top <- merged_signi_drug_summary[which(indiv_drugs_score$medthod_freq >= 2),]
rownames(indiv_drugs_res_top) <- 1:nrow(indiv_drugs_res_top)

write_tsv(indiv_drugs_res_top,file = paste0(output_dir,prefix,"_indiv_top_drugs.tsv"))

