# function to summarize pathways across disease signatures
## add functions for pathway clusters and selecting the top terms based on variance

# last modified: 10/14/25
# Kewalin Samart

library(rrvgo)
library(pheatmap)
library(readr)
library(tidyverse)
library(simplifyEnrichment)

get_GO_asinh_matrix <- function(data_path, pattern, prefix_sub, suffix_sub){
  filenames <- list.files(path = data_path, pattern = pattern, full.names = FALSE)

  GO_mat <- NULL  # initialize result matrix

  for(i in seq_along(filenames)){
    file <- filenames[i]
    GO_res <- read_tsv(file.path(data_path, file), show_col_types = FALSE)
    #GO_res <- GO_res[GO_res$qvalue < 0.05,] # we are filtering out non-significant terms later

    # get signature name
    signature_name <- gsub(paste0(prefix_sub, "_"), "", file)
    signature_name <- gsub(paste0("_", suffix_sub), "", signature_name)
    signature_name <- gsub(".tsv", "", signature_name)
    signature_name <- gsub("_limma", "", signature_name)
    signature_name <- gsub("_DESeq2", "", signature_name)
    signature_name <- make.unique(signature_name)  # ensure unique signature names

    if(nrow(GO_res) > 0){
      # compute LogRatio numeric
      GO_res$asinh_numeric <- as.numeric(GO_res$asinh_observed_expected)

      GO_col <- GO_res[c("Description", "asinh_numeric")]
      colnames(GO_col)[2] <- signature_name

    } else {
      # handle empty result: create single NA row
      GO_col <- data.frame(Description = NA, temp_col = NA)
      colnames(GO_col)[2] <- signature_name
    }

    # merge or initialize GO_mat
    if(is.null(GO_mat)){
      GO_mat <- GO_col
    } else {
      GO_mat <- merge(GO_mat, GO_col, by = "Description", all = TRUE)
    }
  }

  # fill missing values with 0.0 or NA
  if(!is.null(GO_mat)){
    GO_mat[is.na(GO_mat)] <- NA
  }
  GO_mat <- GO_mat[!is.na(GO_mat$Description), ]

  return(GO_mat)
}

get_GO_qvalue_matrix <- function(data_path, pattern, prefix_sub, suffix_sub){
  filenames <- list.files(path = data_path, pattern = pattern, full.names = FALSE)

  GO_mat <- NULL  # initialize result matrix

  for(i in seq_along(filenames)){
    file <- filenames[i]
    GO_res <- read_tsv(file.path(data_path, file), show_col_types = FALSE)

    # get signature name
    signature_name <- gsub(paste0(prefix_sub, "_"), "", file)
    signature_name <- gsub(paste0("_", suffix_sub), "", signature_name)
    signature_name <- gsub(".tsv", "", signature_name)
    signature_name <- gsub("_limma", "", signature_name)
    signature_name <- gsub("_DESeq2", "", signature_name)
    signature_name <- make.unique(signature_name)  # ensure unique signature names

    if(nrow(GO_res) > 0){
      # compute LogRatio numeric
      GO_res$qvalue_numeric <- as.numeric(GO_res$qvalue)

      GO_col <- GO_res[c("Description", "qvalue_numeric")]
      colnames(GO_col)[2] <- signature_name

    } else {
      # handle empty result: create single NA row
      GO_col <- data.frame(Description = NA, temp_col = NA)
      colnames(GO_col)[2] <- signature_name
    }

    # merge or initialize GO_mat
    if(is.null(GO_mat)){
      GO_mat <- GO_col
    } else {
      GO_mat <- merge(GO_mat, GO_col, by = "Description", all = TRUE)
    }
  }

  # fill missing values with 0.0 or NA
  if(!is.null(GO_mat)){
    GO_mat[is.na(GO_mat)] <- NA
  }
  GO_mat <- GO_mat[!is.na(GO_mat$Description), ]

  return(GO_mat)
}

get_GO_testStat_matrix <- function(data_path, pattern, prefix_sub, suffix_sub){
  # data_path <- "/Users/kewalinsamart/Documents/GitHub/drugrep_tb/data/pathways/RNAseq/dn/GO_ORA"
  # prefix_sub <- "GO_ORA_BGcorrected"
  # suffix_sub <- "none.tsv"

  filenames <- list.files(path = data_path,pattern = pattern, all.files = FALSE,
                          full.names = FALSE, recursive = FALSE,ignore.case = FALSE,
                          include.dirs = FALSE, no.. = FALSE)

  for(i in 1:length(filenames)){
    file = filenames[i]
    GO_res <- read_tsv(paste0(data_path,"/",file),show_col_types = FALSE)

    if(nrow(GO_res) > 0){
      # get only data description i.e., signature name
      signature_name <- gsub(paste0(prefix_sub,"_"),"",file)
      signature_name <- gsub(paste0("_",suffix_sub),"",signature_name)

      GO_res$asinh_numeric <- as.numeric(GO_res$asinh_observed_expected)
      GO_col <- GO_res[c("Description","asinh_numeric")]
      colnames(GO_col)[2] <- signature_name

      if(!exists("GO_mat")){
        #union_GO_res <- GO_res
        GO_mat <- GO_col
      }else{
        #union_GO_res <- rbind(union_GO_res, GO_res)
        GO_mat <- merge(GO_mat,GO_col, by = "Description",all = TRUE)
      }
    }
  }
  GO_mat[is.na(GO_mat)] <- 0.0

  return(GO_mat)
}

get_enriched_GOmat <- function(GO_mat, enriched_pct=0.5){
  # this step will remove some terms that are important but not enriched across signatures
  # we can skip this step, and instead pick the lowest adjusted p-value across all enriched GO result tables
  enriched_GO_mat <- GO_mat[rowSums(GO_mat == 0) <= (1-enriched_pct)*ncol(GO_mat),]
  return(enriched_GO_mat)
}

get_GO_annotation <- function(annotation_file='./data/metadata/direct-annotations__Homo_sapiens__Entrez.tsv'){
  GO_annotation = read_tsv(annotation_file,show_col_types = FALSE)
  return(GO_annotation)
}

get_GOterms_cl <- function(pathway_variance_df){
  # map the GO descriptions to GO ids
  GO_df <- merge(pathway_variance_df, GO_annotation, by.x="terms", by.y="Term", all.x=T, all.y=F)
  GO_clean <- na.omit(GO_df)
  GOids_terms <- unique(GO_clean[,c("GOID","terms")])
  GOids <- GOids_terms$GOID
  # calculate similarity matrix
  simMat <- GO_similarity(go_id = GOids, ont = 'BP', db = "org.Hs.eg.db",
                          measure = "Rel",remove_orphan_terms = FALSE)
  # for each GO term, select the highest variance
  cl = cluster_terms(simMat)
  GOids_terms$cluster = cl
  cl_GOids_terms = merge(GOids_terms, pathway_variance_df,by="terms")
  cl_GOids_terms = cl_GOids_terms[order(cl_GOids_terms$cluster, decreasing = FALSE),]

  return(cl_GOids_terms)
}

get_topN_GOterms_cl <- function(cl_GOids_terms, top_by = "var", N = 1){
  if(top_by == "var"){
    top_terms_per_cluster <- cl_GOids_terms %>%
      group_by(cluster) %>%
      slice_max(order_by = variance, n = N, with_ties = FALSE) %>%
      ungroup()
  }else if(top_by == "con"){
    top_terms_per_cluster <- cl_GOids_terms %>%
      group_by(cluster) %>%
      slice_max(order_by = consistent_count, n = N, with_ties = FALSE) %>%
      ungroup()
  }

  return(top_terms_per_cluster)
}
