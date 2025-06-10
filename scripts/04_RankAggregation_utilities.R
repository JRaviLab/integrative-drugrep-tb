# utility functions for #drugrep_tb project
# created date: 01/02/24
# last modified: 09/14/24
# Kewalin Samart

get_full_ranked_drugs <- function(score, drugscore_matrix, ranked_by="median"){
  #' @description this function summarize drugscore_matrix into a full-ranked drug list for a given score sorted by a statistics value of choice
  #' @param score a string indicating one of the connectivity scores: CMAP, WCS, NCS, Tau, Cor_pearson, Cor_spearman
  #' @param drugscore_matrix a drug score matrix in data frame format with first column: unique_pert (drug names), other column names: disease signatures, entries filled with disease-drug scores
  #' @param ranked_by a string indicating the statistics to be used for ranking the drugs
  #' @returns score_ranked_drugs a dataframe of ranked drug names by the chosen ranked_by parameter; column name: score
  #' @author Kewalin Samart

  if(ranked_by=="median"){
    drugscore_matrix_wstats <- drugscore_matrix %>% rowwise() %>%
      mutate(median = median(c_across(where(is.numeric)), na.rm = TRUE)) %>%
      arrange(median)
  }else if(ranked_by=="mean"){
    drugscore_matrix_wstats <- drugscore_matrix %>% rowwise() %>%
      mutate(mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
      arrange(mean)
  }

  score_ranked_drugs <- as.data.frame(drugscore_matrix_wstats[1]) # grab ranked drug names
  colnames(score_ranked_drugs)[1] <- score

  return(score_ranked_drugs)
}

get_partial_ranked_druglists <- function(full_ranked_drugs_matrix, top_pct_drugs_list){
  #' @description this function filter full_ranked_drugs_matrix by the drugs by given sets of top drugs
  #' @param full_ranked_drugs_matrix a dataframe
  #' @param top_pct_drugs_list list of character vectors
  #' @returns partial_ranked_drug_lists a list of vectors of ranked-top drugs by score
  #' @author Kewalin Samart

  partial_ranked_drug_list <- list()
  scores <- colnames(full_ranked_drugs_matrix)
  for(score in scores){
    top_pct_drugs <- top_pct_drugs_list[[score]]
    partial_ranked_drugs <- full_ranked_drugs_matrix[score][full_ranked_drugs_matrix[score][[score]] %in% top_pct_drugs,]
    if(length(partial_ranked_drugs) > 0){
      partial_ranked_drug_list[[score]] <- partial_ranked_drugs
    }
  }
  return(partial_ranked_drug_list)
}

matrix_transfer <- function (ranked_drug_list, full = FALSE){
  #' @description this function transform ranked drug list to a rank matrix, the input format for Big_diffuse
  #' @param druglist ranked list of drugs (list of ranked character vectors filled with drug names) to convert to matrix format
  #' @param full a boolean indicating whether to aggregate the full list or just partial (allows nonoverlapping elements)
  #' @returns rank_mat a matrix containing rank of each element (drug names on the row names)
  #' @source https://github.com/baillielab/comparison_of_RA_methods
  #' @author Kewalin Samart (contribution: modified input arguments)

  unique_nammes = unique(c(ranked_drug_list, recursive = TRUE))

  rank_matrix = matrix(NA, nrow = length(unique_nammes), ncol = length(ranked_drug_list),
                dimnames = list(unique_nammes, names(ranked_drug_list)))

  if(!full){ # partial aggregation: assign the last rank to nonoverlapping elements
    for (i in 1:length(ranked_drug_list)) {
      rank_matrix[, i] = (1 + length(ranked_drug_list[[i]]))
    }
  }
  for (i in 1:length(ranked_drug_list)) {
    rank_matrix[ranked_drug_list[[i]], i] = (1:length(ranked_drug_list[[i]]))
  }
  return(rank_matrix)
}

apply_BiG_RankAggregation <- function(ranked_drug_list, full=FALSE, n_p1=0, M=2000, burnin=1000, prior="IG"){
  #' @description this function applies BiG method implemented with diffuse Inverse Gamma prior or diffuse Uniform prior
  #' for the variance/standard deviation parameters to perform rank aggregation on the input ranked drug list (two platform per aggregation)
  #' @param ranked_drug_list a list of character vector with drug names sorted by statistic of choice
  #' @param full a boolean indicating whether to aggregate the full lists (TRUE) or just partial (FALSE)
  #' @param n_p1 number of studies belong to platform 1
  #' @param M number of MCMC iterations
  #' @param burnin number of burn-in iterations
  #' @param prior BiG mthod arguments
  #' @returns RA_result_df a dataframe with drug names (drug_name) and their rank scores (rank_score)
  #' @source https://github.com/baillielab/comparison_of_RA_methods/blob/main/algorithms/useBiG.R
  #' @author Kewalin Samart (contribution: combining results and entities into data frame)

  # transform to a rank matrix
  rank_matrix <- matrix_transfer(ranked_drug_list, full)

  NTlength <- rep(1000,length(rank_matrix ))
  for (i in 1:length(rank_matrix )) {
    NTlength[i] <- length(rank_matrix[[i]])
  }
  result <- BiG_diffuse(r=rank_matrix, n_T=NTlength, n_p1, M, burnin, prior)
  entities <- rownames(rank_matrix)
  rankedEntities <- entities[order(result,decreasing=TRUE)]

  RA_result_df <- data.frame("drug_name"=entities,"rank_score"=result)
  RA_result_df <- RA_result_df[order(-RA_result_df$rank_score),]

  return(RA_result_df)
}
