# script to perform rank aggregation on the ranked drug results
# Used approach: BiG (A Bayesian latent variable approach i.e., Bayesian Aggregation in Genomic applications)
# **for partial and top ranked lists; incorporating the effect of clustering**
# ref1: https://doi.org/10.1002%2Fsim.7920
# ref2: https://github.com/baillielab/comparison_of_RA_methods/blob/main/algorithms/useBiG.R
# created date: 01/02/24
# last modified: 09/14/24
# Kewalin Samart

source("./scripts/04_RankAggregation_utilities.R")
source("./scripts/04_BiG_RankAggregation_functions.R")

library(tidyverse)
library(readr)

data_technologies <- c("microarray","RNAseq")
n_p1_params <- c(4,3) # number of studies (scores) belong to platform 1 (GSEA-based); RNAseq is missing WCS

## Individual signatures
for(i in range(1:2)){
  data_technology <- data_technologies[i]
  n_p1 <- n_p1_params[i]

  dirname <- paste0("./results/uniformly_processed/",data_technology,"/")
  # read in ranked drug list
  ranked_drug_list <- readRDS(paste0(dirname,"04_rank_aggregation/ranked_0.9pct0.5rv_topdrugs_indivSig_TB_",data_technology,"_list.rds"))

  # transform to a rank matrix
  rank_matrix <- matrix_transfer(ranked_drug_list, full = FALSE)

  NTlength <- rep(1000,length(rank_matrix ))
  for (i in 1:length(rank_matrix )) {
    NTlength[i] <- length(rank_matrix[[i]])
  }

  # perform BiG aggregation approach
  result <- BiG_diffuse(r=rank_matrix, n_T=NTlength, n_p1, M=2000, burnin=1000, prior="IG")
  entities <- rownames(rank_matrix) # obtain drug names
  rankedEntities <- entities[order(result,decreasing=TRUE)] # sort drug names by their rank scores
  RA_0.9pct0.5rv <- data.frame("drug_name"=entities,"rank_score"=result)
  RA_0.9pct0.5rv <- RA_0.9pct0.5rv[order(-RA_0.9pct0.5rv$rank_score),] # put higher rank scores at the top

  # save RA outputs
  write_tsv(RA_0.9pct0.5rv,paste0(dirname,"04_rank_aggregation/RAresult_0.9pct0.5rv_indivSig_TB_",data_technology,".tsv")) # data frame
  saveRDS(rankedEntities, paste0(dirname,"04_rank_aggregation/RAresult_0.9pct0.5rv_indivSig_TB_",data_technology,".rds")) # vector of ranked drugs

  # filter by top methods-wide drugs (picked up by 2/3 method categories: CMAP, LINCS, and Correlation-based methods)
  methodswide_drugs <- read.delim(paste0(dirname,"03b_uni_",data_technology,"_indiv_top_drugs.tsv"),sep = "\t")
  RA_methodswide_drugs <- methodswide_drugs[methodswide_drugs$significant_drug %in% rankedEntities,]
  # save RA methods-wide drug data frame
  write_tsv(RA_methodswide_drugs,paste0(dirname,"04_rank_aggregation/04_uni_",data_technology,"_indiv_RAmethodswide_drugs.tsv"))
}

n_p1_params <- c(4,4)
## Aggregated signatures
for(i in range(1:2)){
  data_technology <- data_technologies[i]
  n_p1 <- n_p1_params[i]

  dirname <- paste0("./results/uniformly_processed/",data_technology,"/")
  # read in ranked drug list
  ranked_drug_list <- readRDS(paste0(dirname,"04_rank_aggregation/ranked_aggrSig_drugs_scores_neg_left_0.9_",data_technology,"_list.rds"))

  # transform to a rank matrix
  rank_matrix <- matrix_transfer(ranked_drug_list, full = TRUE)
  # correct ranks
  rank_matrix_adjusted = as.data.frame(lapply(as.data.frame(rank_matrix), function(x) rank(x, ties.method = "min")))
  row.names(rank_matrix_adjusted) = row.names(rank_matrix)
  rank_matrix_adjusted = as.matrix(rank_matrix_adjusted)

  NTlength <- rep(1000,length(rank_matrix_adjusted))
  for (i in 1:length(rank_matrix_adjusted)) {
    NTlength[i] <- length(rank_matrix_adjusted[[i]])
  }

  # perform BiG aggregation approach
  result <- BiG_diffuse(r=rank_matrix_adjusted, n_T=NTlength, n_p1, M=2000, burnin=1000, prior="IG")
  entities <- rownames(rank_matrix) # obtain drug names
  rankedEntities <- entities[order(result,decreasing=TRUE)] # sort drug names by their rank scores
  RA_negleft0.9 <- data.frame("drug_name"=entities,"rank_score"=result)
  RA_negleft0.9 <- RA_negleft0.9[order(-RA_negleft0.9$rank_score),] # put higher rank scores at the top

  # save RA outputs
  write_tsv(RA_negleft0.9,paste0(dirname,"04_rank_aggregation/RAresult_negleft0.9_aggrSig_TB_",data_technology,".tsv")) # data frame
  saveRDS(rankedEntities, paste0(dirname,"04_rank_aggregation/RAresult_negleft0.9_aggrSig_TB_",data_technology,".rds")) # vector of ranked drugs

  # filter by top methods-wide drugs (picked up by 2/3 method categories: CMAP, LINCS, and Correlation-based methods)
  methodswide_drugs <- read.delim(paste0(dirname,"03b_uni_",data_technology,"_aggr_top_drugs.tsv"),sep = "\t")
  RA_methodswide_drugs <- methodswide_drugs[methodswide_drugs$significant_drug %in% rankedEntities,]
  # save RA methods-wide drug data frame
  write_tsv(RA_methodswide_drugs,paste0(dirname,"04_rank_aggregation/04b_uni_",data_technology,"_aggr_RAmethodswide_drugs.tsv"))
}


