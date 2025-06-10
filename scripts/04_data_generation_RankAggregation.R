# script to generate intermediate data objects for rank aggregation (for results from both indiv and aggr signatures)
# created date: 01/02/24
# last modified: 09/14/24
# Kewalin Samart

library(tidyverse)
library(readr)

# Individual signatures
# generate full ranked drugs matrix -- individual signatures
scores <- c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")
# compute mean and median of each drug
for(data_technology in c("microarray","RNAseq")){
  dirname <- paste0("./results/uniformly_processed/",data_technology,"/")
  i = 1
  for(score in scores){
    drugs_score_mean_med_added <- readRDS(paste0(dirname,score,"_indivSig_TB_",data_technology,"_matrix.rds"))
    score_ranked_drugs <- as.data.frame(drugs_score_mean_med_added$unique_pert)
    colnames(score_ranked_drugs)[1] <- score
    if(i > 1){
      score_ranked_drugs_df <- cbind(score_ranked_drugs_df,score_ranked_drugs)
    }else{
      score_ranked_drugs_df <- score_ranked_drugs
    }
    i = i + 1
  }
  saveRDS(score_ranked_drugs_df, paste0(dirname,"04_rank_aggregation/full_ranked_drugs_indivSig_TB_",data_technology,"_matrix.rds"))
}

# generate lists of partial ranked top drugs and unranked top drugs from pct approach
top_pct_drugs_list <- list()
score_ranked_drugs_list <- list()
scores <- c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")
for(data_technology in c("microarray","RNAseq")){
  dirname <- paste0("./results/uniformly_processed/",data_technology,"/")
  for(score in scores){
    top_drugs_tech_df <- read.delim(paste0("./results/uniformly_processed/",data_technology,"/",score,"_top_0.9pct_0.5reversed_drugs.tsv"),sep = "\t")
    drugs_score_mean_med_added <- readRDS(paste0(dirname,score,"_indivSig_TB_",data_technology,"_matrix.rds")) # overwrite the current matrix with same data and mean/median added
    topdrugs_score_mean_med_added <- drugs_score_mean_med_added[drugs_score_mean_med_added$unique_pert %in% top_drugs_tech_df$unique_pert,]
    score_ranked_drugs <- topdrugs_score_mean_med_added$unique_pert
    top_pct_drugs_list[[score]] <- top_drugs_tech_df$unique_pert
    score_ranked_drugs_list[[score]] <- score_ranked_drugs
  }
  saveRDS(score_ranked_drugs_list, paste0(dirname,"04_rank_aggregation/ranked_0.9pct0.5rv_topdrugs_indivSig_TB_",data_technology,"_list.rds"))
  saveRDS(top_pct_drugs_list, paste0(dirname,"04_rank_aggregation/unranked_0.9pct0.5rv_topdrugs_indivSig_TB_",data_technology,"_list.rds"))
}

#-----------------
# Aggregated signatures
# generate lists of partial ranked top drugs
for(data_technology in c("microarray","RNAseq")){
  dirname <- paste0("./results/uniformly_processed/",data_technology,"/")
  drugs_score <- read.csv(paste0(dirname,"uni_",data_technology,"_aggrSig_drugs_scores_neg_left_0.9.tsv"),sep="\t")
  for(i in 3:dim(drugs_score)[2]){
    score_ranked_drugs = drugs_score[order(drugs_score[i], decreasing = FALSE), ]$pert
    score_ranked_drugs <- as.data.frame(score_ranked_drugs)
    if(i > 3){
      score_ranked_drugs_df <- cbind(score_ranked_drugs_df,score_ranked_drugs)
    }else{
      score_ranked_drugs_df <- score_ranked_drugs
    }
  }
  colnames(score_ranked_drugs_df) <- c("CMAP","Tau","NCS","WCS","Cor_spearman","Cor_pearson")
  score_ranked_drugs_df <- score_ranked_drugs_df[c("CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson")]
  score_ranked_drugs_list <- lapply(names(score_ranked_drugs_df), function(col) score_ranked_drugs_df[[col]])
  names(score_ranked_drugs_list) <- colnames(score_ranked_drugs_df)
  saveRDS(score_ranked_drugs_df, paste0(dirname,"04_rank_aggregation/full_ranked_drugs_aggrSig_TB_",data_technology,"_matrix.rds"))
  saveRDS(score_ranked_drugs_list, paste0(dirname,"04_rank_aggregation/ranked_aggrSig_drugs_scores_neg_left_0.9_",data_technology,"_list.rds"))
}

