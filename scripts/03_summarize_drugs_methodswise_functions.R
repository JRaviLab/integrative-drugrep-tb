# functions to finalize drug results from individual disease signatures
# created dates: 06/13/22
# last modified: 09/08/24
# Kewalin Samart

clean_data_to_run <- function(metadata_path,disease_keywords){
  #' @description this function subsets metadata of disease signatures of interest i.e., selecting specific disease comparisons
  #' @param metadata_path path to metadata file of input disease signatures
  #' @param disease_keywords a string indicating specific disease keywords separated by comma
  #' @returns data_to_run: filtered metadata of disease signatures of interest
  #' @author Kewalin Samart

  # read in metadata file
  data_to_run <- read.delim(metadata_path, sep = "\t")
  # subset the metadata table to get only data for the specified disease
  # subset specific infection(s)
  data_to_run <- data_to_run[grep(paste(disease_keywords, collapse = "|"),
                                  data_to_run$file_name), ]
  # subset data table to get only infected-control signatures
  data_to_run <- data_to_run[grep("control", data_to_run$file_name), ]

  return(data_to_run)
}

get_drug_results <- function(data_to_run, drug_res_path, score_method, score){
  #' @description function to read in drug results prioritized by a specific connectivity score
  #' @param data_to_run metadata for reading in drug results containing the columns: "accession_no" indicating GEO/refine.bio accession id e.g., GSE16250
  #' @param data_to_run (cont.) and "file_name"indicating cond1_cond2 (DE comparison for the disease signature) e.g., MTB_control
  #' @param drug_res_path path to the data folder with drug results
  #' @param score_method a string indicating a choice of connectivity methods: "CMAP","LINCS","Cor_spearman","Cor_pearson"
  #' @param score a string indicating a choice of connectivity scores: "CMAP","WCS","NCS","Tau","Cor_spearman","Cor_pearson"
  #' @returns combined_drug_df: a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @author Kewalin Samart
  require(stringr)

  for(i in 1:nrow(data_to_run)){
    # read in variables
    accession_no <- data_to_run$accession_no[i] # GEO/refine.bio accession id e.g., GSE16250
    filename <- data_to_run$file_name[i] # cond1_cond2 e.g., MTB_control

    #print(paste0("Iteration ",i))
    #print(paste0("Getting drug result: ", accession_no," ", filename))

    drug_df <- read.delim(paste0(drug_res_path,score_method,"_",accession_no,"_",filename,".tsv"), sep="\t")
    #drug_df <- merge(x=drug_df, y=drug_metadata, by.x = "pert", by.y = "pert_iname", all.x=TRUE, all.y=FALSE) # not needed as we are merging drug info later

    if(score == "CMAP"){
      drug_df <- drug_df[, c("scaled_score","pert","cell","t_gn_sym")]
    }else if(score == "WCS"){
      drug_df <- drug_df[, c("WTCS","pert","cell","t_gn_sym")]
    }else if(score == "NCS"){
      drug_df <- drug_df[, c("NCS","pert","cell","t_gn_sym")]
    }else if(score == "Tau"){
      drug_df <- drug_df[, c("Tau","pert","cell","t_gn_sym")]
    }else if(score %in% list("Cor_spearman","Cor_pearson")){
      drug_df <- drug_df[, c("cor_score","pert","cell","t_gn_sym")]
    }
    # add pert_cell column
    drug_df$pert_cell <- str_c(drug_df$pert,"_",drug_df$cell)

    # add GSE_platform_contrast column
    drug_df$GSE_platform <- paste0(accession_no,"_",filename)

    if(i == 1){
      combined_drug_df <- drug_df
    }else{
      # bind dataframes by rows -- rbind
      combined_drug_df <- rbind(drug_df, combined_drug_df)
    }
  }
  return(combined_drug_df)
}

get_DrugDis_ScoreMatrix <- function(combined_drug_df, score, stats = "median"){
  #' @description this function converts the combined_drug_df to a dataframe (matrix filled with scores) where rows: drugs and cols: disease signatures
  #' @param combined_drug_df a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @param score a string indicating a string indicating a choice of connectivity scores: "CMAP","WCS","NCS","Tau","Cor_spearman", "Cor_pearson"
  #' @param stats a string indicating the stats of score to include in the matrix: "median", "min", max".; "median" by default.
  #' @returns transformed_combined_drug_df: a dataframe with rows: drugs, columns: disease signatures, entries: scores
  #' @author Kewalin Samart

  # modify the dataframe for filtering significant drugs
  unique_pert <- unique(combined_drug_df$pert)
  unique_signature <- unique(combined_drug_df$GSE_platform)

  transformed_combined_drug_df <- as.data.frame(unique_pert)

  col_idx = 1
  for(sig in unique_signature){
    col_idx = col_idx + 1
    # get tau vector for each signature
    score_vec = c()
    for(pert in unique_pert){

      pert_df = combined_drug_df[which(combined_drug_df$pert == pert & combined_drug_df$GSE_platform == sig),]
      if(score == "Tau"){
        if(stats == "min"){
          index = which.min(pert_df$Tau)
          score_vec <- append(score_vec, pert_df$Tau[index])
        }else if(stats == "median"){
          score_vec <- append(score_vec, median(pert_df$Tau))
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$Tau))
        }
      }
      if(score == "NCS"){
        if(stats == "min"){
          index = which.min(pert_df$NCS)
          score_vec <- append(score_vec, pert_df$NCS[index])
        }else if(stats == "median"){
          score_vec <- append(score_vec, median(pert_df$NCS))
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$NCS))
        }
      }
      if(score == "WCS"){
        if(stats == "min"){
          index = which.min(pert_df$WTCS)
          score_vec <- append(score_vec, pert_df$WTCS[index])
        }else if(stats == "median"){
          score_vec <- append(score_vec, median(pert_df$WTCS))
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$WTCS))
        }
      }
      if(score == "CMAP"){
        if(stats == "min"){
          index = which.min(pert_df$scaled_score)
          score_vec <- append(score_vec, pert_df$scaled_score[index])
        }else if(stats == "median"){
          score_vec <- append(score_vec, median(pert_df$scaled_score))
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$scaled_score))
        }
      }
      if(score == "Cor_spearman" | score == "Cor_pearson"){
        if(stats == "min"){
          index = which.min(pert_df$cor_score)
          score_vec <- append(score_vec, pert_df$cor_score[index])
        }else if(stats == "median"){
          score_vec <- append(score_vec, median(pert_df$cor_score))
        }else if(stats == "max"){
          score_vec <- append(score_vec, max(pert_df$cor_score))
        }
      }
    }
    # now we have the vector; we want to add this as a new column to the data frame
    transformed_combined_drug_df[col_idx] <- score_vec
    colnames(transformed_combined_drug_df)[col_idx] <- sig
  }

  return(transformed_combined_drug_df)
}

get_DrugCellDis_ScoreMatrix <- function(combined_drug_df, score){
  #' @description this function converts the combined_drug_df to a dataframe (matrix filled with scores) where rows: drug-cell combinations and cols: disease signatures
  #' @param combined_drug_df a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @param score a string indicating a string indicating a choice of connectivity scores: "CMAP","WCS","NCS","Tau","Cor_spearman", "Cor_pearson"
  #' @returns transformed_combined_drug_df: a dataframe with rows: drug-cell combinations, columns: disease signatures, entries: scores
  #' @author Kewalin Samart

  # modify the dataframe for filtering significant drugs; considering cell types
  unique_pert_cell <- unique(combined_drug_df$pert_cell)
  unique_signature <- unique(combined_drug_df$GSE_platform)

  transformed_combined_drugcell_mat = matrix(NA, nrow = length(unique_pert_cell), ncol = length(unique_signature),
                                             dimnames = list(unique_pert_cell, unique_signature))
  col_idx = 1
  for(sig in unique_signature){
    row_idx = 1
    for(pert_cell in unique_pert_cell){
      pert_cell_df = combined_drug_df[which(combined_drug_df$pert_cell == pert_cell & combined_drug_df$GSE_platform == sig),]
      if(nrow(pert_cell_df) != 0){
        if(score == "Tau"){
          score_val = pert_cell_df$Tau[1]
        }
        if(score == "NCS"){
          score_val = pert_cell_df$NCS[1]
        }
        if(score == "WCS"){
          score_val = pert_cell_df$WTCS[1]
        }
        if(score == "CMAP"){
          score_val = pert_cell_df$scaled_score[1]
        }
        if(score == "Cor_spearman" | score == "Cor_pearson"){
          score_val = pert_cell_df$cor_score[1]
        }
        transformed_combined_drugcell_mat[row_idx,col_idx] = score_val
      }else{
        transformed_combined_drugcell_mat[row_idx,col_idx] = 0
      }
      row_idx = row_idx + 1
    }
    col_idx = col_idx + 1
  }
  transformed_combined_drugcell_df <- as.data.frame(transformed_combined_drugcell_mat)
  transformed_combined_drugcell_df$unique_pert <- row.names(transformed_combined_drugcell_df)
  transformed_combined_drugcell_df <- transformed_combined_drugcell_df[c("unique_pert",unique_signature)]
  row.names(transformed_combined_drugcell_df) <- NULL

  return(transformed_combined_drugcell_df)
}

determine_threshold <- function(score_vec, values, tail,  score_percentile=0.90){
  #' @description this function computes a threshold for getting top drugs by a predefined score percentile
  #' @param score_vec a numeric vector of disease-drug scores
  #' @param score_percentile a numeric indicating percentile of choice e.g. 0.75 representing a threshold being the 75th percentile of the flipped score distribution (dist x-axis: positive <---> negative, as negative score preferred for disease-drug reversal)
  #' @param values a string indicating which non-zero values to include in the distribution: "neg", "pos", if other values specified for this argument, it means including all values.
  #' @param tail a string indicating which side of the distribution to get the indicated top percentile from
  #' @returns threshold
  #' @author Kewalin Samart

  # set scores of zero to NA as they are meaningless
  score_vec[score_vec == 0] <- NA
  # considerations of what non-zero values to include in threshold calculation
  if(values == "neg"){
    print("Only consider negative values")
    score_vec[score_vec > 0] <- NA # set positive values to NAs so they would be removed later
  }else if(values == "pos"){
    print("Only consider positive values")
    score_vec[score_vec < 0] <- NA # set positive values to NAs so they would be removed later
  }

  score_vec <- score_vec[!is.na(score_vec)] # remove NAs

  # compute threshold based on a given percentile and tail direction
  if(tail == "left"){
    threshold <- round(unname(quantile(score_vec, 1-score_percentile)[paste0((1-score_percentile)*100,"%")])[1],3)
  }else if(tail == "right"){
    threshold <- round(unname(quantile(score_vec, score_percentile)[paste0((score_percentile)*100,"%")])[1],3)
  }

  print(paste0("threshold:",threshold))

  return(threshold)
}

get_reversing_drugs_freq <- function(transformed_combined_drug_df, threshold, values){
  #' @description this function counts how often the drug got prioritized by the input disease signatures based on a given score threshold
  #' @param transformed_combined_drug_df a dataframe with rows: drugs, columns: disease signatures, entries: scores
  #' @param threshold a negative numeric indicating the threshold for reversal
  #' @param values a string indicating which non-zero values to include in the distribution: "neg", "pos", if other values specified for this argument, it means including all values.
  #' @returns pert_occurrence_df: a sorted dataframe with drug names and their occurrences across the input signatures
  #' @author Kewalin Samart

  # get unique drug names and initialize a dataframe for output
  unique_pert <- unique(transformed_combined_drug_df$unique_pert)
  pert_occurrence_df <- as.data.frame(unique_pert)

  # prepare a matrix-like dataframe for counting occurrences
  rownames(transformed_combined_drug_df) <- transformed_combined_drug_df$unique_pert
  transformed_combined_drug_df$unique_pert <- NULL

  # loop through each drug and disease signatures to get occurrence counts
  pert_above_occurrence <- list()
  for(row in 1:nrow(transformed_combined_drug_df)){
    above_count <- 0
    for(col in 1:ncol(transformed_combined_drug_df)){
      if(values == 'neg'){
        if(transformed_combined_drug_df[row,col] < threshold){
          above_count = above_count + 1
        }
      }else if(values == 'pos'){
        if(transformed_combined_drug_df[row,col] > threshold){
          above_count = above_count + 1
        }
      }else if(values == 'zero'){
        if(abs(transformed_combined_drug_df[row,col]) < threshold){
          above_count = above_count + 1
        }
      }
    }
    pert_above_occurrence <- append(pert_above_occurrence, above_count)
  }
  # add occurrence column to the dataframe
  pert_occurrence_df$occurrence <- as.numeric(pert_above_occurrence)

  # add drug info to the occurrence dataframe
  # "moa","disease_area","t_gn_sym","target","indication"
  # read in drug info from drug repurposing hub (treat this as a extdata)
  drug_metadata <- read.csv(file="./data/metadata/repurposing_drugs_20200324.csv",skip = 9)
  pert_occurrence_df <- merge(x=pert_occurrence_df, y=drug_metadata, by.x = "unique_pert", by.y = "pert_iname", all.x=TRUE, all.y=FALSE)
  # select drugs with at least 1 occurrence
  pert_occurrence_df <- pert_occurrence_df[pert_occurrence_df$occurrence > 0,]
  # sort drugs by their occurrence in decreasing order
  pert_occurrence_df = pert_occurrence_df[order(pert_occurrence_df$occurrence, decreasing = TRUE),]

  return(pert_occurrence_df)
}


get_significant_drugs <- function(pert_occurrence_df, transformed_combined_drug_df, percent_reverse = NA, n = NA){
  #' @description this function identifies significant drugs based on reversing threshold and either a given percentage of reversal or top n drugs ranked by the score
  #' @param pert_occurrence_df a sorted dataframe with drug names and their occurrences across the input signatures
  #' @param transformed_combined_drug_df a dataframe with rows: drugs, columns: disease signatures, entries: scores
  #' @param percent_reverse a numerical value in the range [0,1] indicating the percentage of reversal for drug prioritization
  #' @param n an integer indicating top n drugs decreasingly sorted by the occurrences
  #' @returns signi_drug_df: a dataframe containing the drugs prioritized "significant" based on either the given percent_reverse or (top) n
  #' @author Kewalin Samart

  # check if either percent_reverse or n is specified to a numeric value
  if(!(is.numeric(percent_reverse) | is.numeric(n))){
    print("Please enter a numerical value for either percent_reverse or n")
    return(NULL)
  }

  # if n (top n drugs) not specified, then get significant drugs by percent of reversal
  if(is.na(n)){
    print("Selecting drugs by percent_reverse")
    signi_drug_df <- pert_occurrence_df[which(pert_occurrence_df$occurrence >= percent_reverse*ncol(transformed_combined_drug_df)),]
  }
  # if percent of reversal not specified, then get significant drugs by top n with the highest occurrences
  if(is.na(percent_reverse)){
    print("Selecting drugs by top n")
    # rank by occurrence frequency
    pert_occurrence_df = pert_occurrence_df[order(pert_occurrence_df$occurrence, decreasing = TRUE),]
    # take the top n drugs if n < number of rows
    if(n < dim(pert_occurrence_df)[1]){
      signi_drug_df <- pert_occurrence_df[1:n,]
    }else{
      signi_drug_df <- pert_occurrence_df
    }

  }
  # check the output dataframe dimension
  if(dim(signi_drug_df)[1] == 0){
    print("No significant drugs identified. Please try adjusting the following parameters: threshold or percent_reverse")
  }
  print(paste0("number of significant drugs:",dim(signi_drug_df)[1]))
  return(signi_drug_df)
}

get_signi_info <- function(combined_drug_df, signi_drug_df){
  #' @description this function grabs the LINCS signatures information of identified significant drugs from the original combined drug dataframe
  #' @param combined_drug_df a dataframe (from multiple disease signatures) with scores, drug names (perturbageons; pert), cell lines of the drug profiles, moas (mechanism of actions), disease areas, and drug target genes
  #' @param signi_drug_df a dataframe containing the drugs prioritized "significant" based on either the given percent_reverse or (top) n
  #' @returns  signi_info_df: a dataframe with significant drugs and their LINCS signatures information
  #' @author Kewalin Samart
  signi_info_df <- combined_drug_df[which(combined_drug_df$pert %in% signi_drug_df$unique_pert),]

  return(signi_info_df)
}

get_signidrugs_by_score <- function(metadata_path, drug_summ_path, dirname="./results"){
  #' @description this function quantifies significant drugs by score given a set of disease signatures in metadata_path
  #' @description (cont.) and score parameters in drug_summ_path: (i) drug_res_path, (ii) score_method, (iii) score, (iv) percent_reverse, (v) n or score_percentile
  #' @param metadata_path path to input disease signatures
  #' @param drug_sum_path path to the drug summarization arguments file
  #' @param dirname output directory; set to "./results/" by default
  #' @returns all_signi_drugs_df: a dataframe containing all significant drugs picked up by one of the 6 scores: CMAP, LINCS: WCS, NCS, Tau, Cor: Cor_spearman, Cor_pearson
  #' @author Kewalin Samart

  #metadata_path <- "./inputs/TB_microarray_args.tsv"
  #drug_summ_path <- "./inputs/drug_summarization_args.tsv"
  #dirname <- paste0("./results/uniformly_processed/microarray/",score_method,"/")

  require(readr)
  # read in metadata for disease signatures and drug results
  data_to_run <- read.delim(file=metadata_path, sep="\t")
  drug_summ_args <- read.delim(file=drug_summ_path, sep="\t")

  all_signi_drugs <- c()
  signi_drugs_score_list <- list()

  for(i in 1:ncol(drug_summ_args)){

    # getting arguments for drug summarization
    drug_res_path <- drug_summ_args$drug_res_path[i]
    score_method <- drug_summ_args$score_method[i]
    score <- drug_summ_args$score[i]
    percent_reverse <- drug_summ_args$percent_reverse[i]
    n <- drug_summ_args$n[i]
    score_percentile <- drug_summ_args$score_percentile[i]

    print(paste("Getting significant drugs for the disease signatures", metadata_path,"prioritized by", score,sep = " "))
    combined_drug_df <- get_drug_results(data_to_run = data_to_run, drug_res_path = drug_res_path, score_method = score_method, score=score)
    threshold <- determine_threshold(combined_drug_df, score_percentile=score_percentile)

    transformed_combined_drug_df <- modify_df_drugs(combined_drug_df,score = score)
    pert_occurrence_df <- get_reversing_drugs_freq(transformed_combined_drug_df, threshold=threshold)
    signi_drug_df <- get_significant_drugs(pert_occurrence_df, transformed_combined_drug_df, percent_reverse = percent_reverse, n = n)

    dirname <- paste(dirname,score_method,sep="/")
    if(!dir.exists(dirname)){
      dir.create(dirname)
    }

    write_tsv(signi_drug_df,file = paste0(dirname,"/",score,"_top_",score_percentile,"pct_",percent_reverse,"reversed_drugs.tsv"))

    signi_drugs_score <- signi_drug_df$unique_pert
    signi_drugs_score_list[[score]] <- signi_drugs_score # named list containing significant drugs prioritized by different metrics
    all_signi_drugs <- c(all_signi_drugs,signi_drugs_score)
  }

  all_signi_drugs_df <- data.frame(unique(all_signi_drugs))
  colnames(all_signi_drugs_df)[1] <- "significant_drug"

  for(score in names(signi_drugs_score_list)){
    score_signi_drugs_boolean <- all_signi_drugs_df$significant_drug %in% signi_drugs_score_list[[score]]
    score_signi_drugs_numeric <- replace(score_signi_drugs_boolean, score_signi_drugs_boolean==TRUE, 1)
    all_signi_drugs_df[score] <- score_signi_drugs_numeric
  }

  return(all_signi_drugs_df)
}

summarize_drugs_bymethods <- function(all_signi_drugs_df, prefix, dirname="./results"){
  #' @description this function prioritizes final significant drugs based on 3 methods: CMAP, LINCS, and Cor
  #' @description (cont.) the drugs that come up in two out of the three methods are selected.
  #' @param all_signi_drugs_df a dataframe containing all significant drugs picked up by one of the 6 scores: CMAP, LINCS: WCS, NCS, Tau, Cor: Cor_spearman, Cor_pearson
  #' @param prefix a string for output file prefix e.g., uni_marray for uniformly processed-microarray data
  #' @param dirname output directory; set to "./results/" by default
  #' @returns indiv_drugs_res_top: a dataframe containing final significant drugs based on 3 methods: CMAP, LINCS, and Cor sorted by their occurrence across the 6 scores
  #' @author Kewalin Samart

  require(readr)
  # summarize based on 3 methods: CMAP, LINCS, Cor
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

  if(!dir.exists(dirname)) {
    dir.create(dirname)
  }

  write_tsv(merged_signi_drug_summary,file = paste0(dirname,"/",prefix,"_indiv_drug_summary.tsv"))

  # replace non-zero values with 1
  indiv_drugs_score <- merged_signi_drug_summary[,c("CMAP","LINCS","Cor")]
  indiv_drugs_score[-1] <- as.integer(indiv_drugs_score[-1] != 0)
  # subset only drugs that show up in at least 2 out of the 3 methods
  indiv_drugs_score$medthod_freq <- rowSums(indiv_drugs_score)
  indiv_drugs_res_top <- merged_signi_drug_summary[which(indiv_drugs_score$medthod_freq >= 2),]
  rownames(indiv_drugs_res_top) <- 1:nrow(indiv_drugs_res_top)

  write_tsv(indiv_drugs_res_top,file = paste0(dirname,"/",prefix,"_indiv_top_drugs.tsv"))

  return(indiv_drugs_res_top)
}

