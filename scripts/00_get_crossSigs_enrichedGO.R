# function to summarize pathways across disease signatures
## change asinh(observed/expected)

# last modified: 10/05/25
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

summarize_unionGOids <- function(GO_mat, GO_annotation){
  union_GOids <- GO_mat$Description
  # map the GO descriptions to GO ids
  union_GO_df <- merge(GO_mat,GO_annotation,by.x="Description",by.y="Term", all.x=T, all.y=F)
  union_GO_clean <- na.omit(union_GO_df)
  union_GOids_terms <- unique(union_GO_clean[,c("GOID","Description")]) # some terms don't have associated GO id
  union_GOids <- union_GOids_terms$GOID
  # calculate similarity matrix
  simMat <- GO_similarity(go_id = union_GOids, ont = 'BP', db = 'org.Hs.eg.db',
                          measure = "Rel",remove_orphan_terms = FALSE)
  # for each GO term, select the lowest adj.p-val
  cl = cluster_terms(simMat)
  union_GOids_terms$cluster = cl
  cl_unionGOids_terms = merge(union_GOids_terms,GO_mat,by="Description")
  cl_unionGOids_terms = cl_unionGOids_terms[order(cl_unionGOids_terms$cluster, decreasing = FALSE),]

  GO_mat_float <- GO_mat[GO_mat$Description %in% union_GOids_terms$Description,]
  GO_mat_float$Description <- NULL
  #min_GOqvals <- apply(GO_mat_float,1,function(val) if(all(val==0)) 0 else min(val[val>0]))
  min_GOqvals <- apply(GO_mat_float, 1, function(val) {
    if (all(val == 0)) 0 else min(val[val > 0])
  })
  cl_unionGOids_terms <- cl_unionGOids_terms %>% mutate(min_qval = min_GOqvals, .after=cluster)
  signatures = colnames(cl_unionGOids_terms)[5:(ncol(cl_unionGOids_terms))]
  cl_unionGOids_terms <- cl_unionGOids_terms[c(c("Description","GOID","min_qval","cluster"),signatures)]

  return(cl_unionGOids_terms)
}

get_clusterGOterms_sigs <- function(cl_unionGOids_terms){
  cluster_num <- unique(cl_unionGOids_terms$cluster)
  signatures <- colnames(cl_unionGOids_terms)[6:ncol(cl_unionGOids_terms)]

  for (cluster in cluster_num) {
    cl_unionGOids_subdf <- cl_unionGOids_terms[cl_unionGOids_terms$cluster == cluster, ]
    sig_cl_qvals <- c()

    for (sig_itr in seq_along(signatures)) {
      #cl_qval <- unname(apply(cl_unionGOids_subdf[5 + sig_itr], 2, function(val) {
      #  if (all(val == 0)) 0 else min(val[val > 0])
      #}))
      cl_qval <- unname(apply(cl_unionGOids_subdf[5 + sig_itr], 2, function(val) {
        val_non_na <- val[!is.na(val)]
        if (length(val_non_na) == 0) {
          return(NA)  # All NA
        } else if (all(val_non_na == 0)) {
          return(0)   # All 0
        } else {
          return(min(val_non_na[val_non_na > 0]))
        }
      }))

      sig_cl_qvals <- c(sig_cl_qvals, cl_qval)
    }

    names(sig_cl_qvals) <- signatures
    sig_cl_qvals_df <- as.data.frame(as.list(sig_cl_qvals))

    if (cluster == cluster_num[1]) {
      combined_df <- sig_cl_qvals_df
    } else {
      combined_df <- rbind(combined_df, sig_cl_qvals_df)
    }
  }

  rownames(combined_df) <- cluster_num
  return(combined_df)
}


get_topthree_sigGOterms_cl <- function(cl_unionGOids_terms){
  # calculate the top three significant GO terms in each cluster
  # get the most significant q-value of each cluster per signature
  cluster_num = unique(cl_unionGOids_terms$cluster)
  signatures = colnames(cl_unionGOids_terms)[4:length(cl_unionGOids_terms)]
  # iterate through each cluster
  for(cluster in cluster_num){
    cl_unionGOids_subdf = cl_unionGOids_terms[cl_unionGOids_terms$cluster == cluster,]
    sorted_by_qval <- cl_unionGOids_subdf[order(cl_unionGOids_subdf$min_qval, decreasing = FALSE),]
    topN_GOterms <- sorted_by_qval$Description[1:3]
    names(topN_GOterms) <- c("top1","top2","top3")
    topN_GOterms_t = t(as.data.frame(topN_GOterms))
    # we need to stack them right here
    if(cluster == 1){
      topN_GOterms_df = topN_GOterms_t
    }else{
      topN_GOterms_df = rbind(topN_GOterms_df, topN_GOterms_t)
    }
  }
  row.names(topN_GOterms_df) = cluster_num
  topN_GOterms_df = as.data.frame(topN_GOterms_df)
  return(topN_GOterms_df)
}



## Example run
# data_path <- "/Users/kewalinsamart/Documents/GitHub/drugrep_tb/data/pathways/microarray/dn/GO_ORA"
# pattern <- "GO_ORA_BGcorrected"
# prefix_sub <- "GO_ORA_BGcorrected"
# suffix_sub <- "none.tsv"
#
#
#
# GO_mat = get_GOqval_matrix(data_path, pattern, prefix_sub, suffix_sub)
# GO_annotation = get_GO_annotation()
# cl_unionGOids_terms = summarize_unionGOids(GO_mat, GO_annotation)
# sig_cl_qvals_mat = get_clusterGOterms_sigs(cl_unionGOids_terms) # get qval matrix: cluster x signature
# write_tsv(sig_cl_qvals_mat,file="./data/pathways/microarray/up/GO_ORA/sig_cl_qvals_mat.tsv")
#
# topN_GOterms_df = get_topthree_sigGOterms_cl(cl_unionGOids_terms) # get top 3 terms enriched in each cluster
# write_tsv(topN_GOterms_df,file="./data/pathways/microarray/up/GO_ORA/topN_GOterms.tsv")
