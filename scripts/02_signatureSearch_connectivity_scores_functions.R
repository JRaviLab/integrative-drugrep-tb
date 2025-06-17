# functions for drug search by CMap1.0, CMap2.0, and Correlation-based connectvity scores
# last modified: 06/17/25
# Kewalin Samart

# loading the needed libraries
library(signatureSearch)
library(rhdf5)
library(ExperimentHub)
library(tidyverse)
library(here)

source(here("scripts/01_signature_landmark_prep_functions.R"))

combine_signatures <- function(up_signature, dn_signature, L1000 = TRUE){
  #' @description This function combines up and dn input signatures into a full signature
  #' @param up_signature up signature dataframe containing genes and other statistics from DE analysis
  #' @param dn_signature dn signature dataframe containing genes and other statistics from DE analysis
  #' @param L1000 a boolean indicating whether to keep only L1000 genes
  #' @returns full_signature
  #' @author Kewalin Samart
  full_signature_ <- rbind(up_signature, dn_signature)
  clean_full_sig = prepare_signature(signature = full_signature_, L1000 = L1000)
  full_signature <- clean_full_sig %>% arrange(-clean_full_sig$log2FoldChange)

  return(full_signature)
}

get_full_signature <- function(full_sig_path, L1000 = TRUE){
  #' @description This function reads in a full signature + convert it to a matrix with row names: GeneID (Entrezid) and a column of log2FoldChange
  #' @param full_sig_path path to a full signature (up+dn genes)
  #' @returns full_signature input for gess_cor methods; signature in the format with row names: GeneID (Entrezid) and a column of log2FoldChange
  #' @param L1000 a boolean indicating whether to keep only L1000 genes
  #' @author Kewalin Samart
  input_full_sig <- read.delim(full_sig_path, sep="\t")
  clean_full_sig = prepare_signature(signature = input_full_sig, L1000 = L1000)
  full_signature <- input_full_sig[,c("GeneID","log2FoldChange")]
  row.names(full_signature) <- as.character(full_signature$GeneID)
  full_signature$GeneID <- NULL
  full_signature <- na.omit(full_signature) # remove gene(s) with NA log2FoldChange
  full_signature_matrix <- as.matrix(full_signature)

  return(full_signature_matrix)
}

get_updn_signature <- function(updn_sig_path, L1000 = TRUE){
  #' @description This function gets up/dn gene lists from either up or down input signatures
  #' @param sig_path path to either up/dn signatures to get genes from
  #' @param L1000 a boolean indicating whether to keep only L1000 genes
  #' @returns updn_genes
  #' @author Kewalin Samart
  input_updn_sig <- read.delim(updn_sig_path, sep='\t')
  clean_input_updn_sig = prepare_signature(signature = input_updn_sig, L1000 = L1000)
  updn_genes <- as.character(clean_input_updn_sig$GeneID)

  return(updn_genes)
}

setup_cmap1db <- function(){
  #' @description  This function sets up CMAP database
  #' @details source code from signatureSearchData package manual
  #' @source  https://bioconductor.org/packages/release/data/experiment/manuals/signatureSearchData/man/signatureSearchData.pdf
  #' @details CMap2 LFC (log2 fold change) Signature Database
  #' @returns cmap_path path to CMAP database
  #' @author Kewalin Samart

  require(rhdf5)
  require(ExperimentHub)
  eh <- ExperimentHub()
  query(eh, c("signatureSearchData", "cmap"))
  cmap_path <- eh[["EH3223"]]
  rhdf5::h5ls(cmap_path)

  return(cmap_path)
}

setup_lincsdb <- function(){
  #' @description This function sets up LINCS database; only contains drug signatures of 10mu and 24h
  #' @details source code from signatureSearchData package manual
  #' @source https://bioconductor.org/packages/release/data/experiment/manuals/signatureSearchData/man/signatureSearchData.pdf
  #' @details LINCS Z-score Signature Database
  #' @returns lincs_path path to LINCS database
  #' @author Kewalin Samart

  require(rhdf5)
  require(ExperimentHub)
  eh <- ExperimentHub()
  query(eh, c("signatureSearchData", "lincs"))
  lincs_path <- eh[["EH3226"]]
  rhdf5::h5ls(lincs_path)

  return(lincs_path)
}

compute_CMap1_scores <- function(up_genes, dn_genes, db_path){
  #' @description This function quantifies drug candidates for a given pair of up/dn ENTREZID vectors prioritized by cmap1
  #' @param up_genes up GeneID (Entrezid; character vector)
  #' @param dn_genes dn GeneID (Entrezid; character vector)
  #' @param db_path obtained form either setup_cmap1db() or setup_lincsdb() functions
  #' @returns final_res a data frame of drug candidates prioritized by cmap1
  #' @author Kewalin Samart

  require(signatureSearch)
  # get CMap 1.0 connectivity scores
  qsig <- qSig(query = list(upset = up_genes,
                            downset = dn_genes),
               gess_method = "CMAP",  # gess_method: one of 'CMAP', 'LINCS', 'gCMAP', 'Fisher' or 'Cor'
               refdb = db_path)
  # obtain query results
  res <- gess_cmap(qsig)
  final_res <- finalize_result(result(res))

  return(final_res)
}

compute_CMap2lincs_scores <- function(up_genes, dn_genes, db_path){
  #' @description This function quantifies drug candidates for a given pair of up/dn ENTREZID vectors prioritized by WCS, NCS, or Tau
  #' @param up_genes up GeneID (Entrezid; character vector)
  #' @param dn_genes dn GeneID (Entrezid; character vector)
  #' @param db_path obtained form either setup_cmap1db() or setup_lincsdb() functions
  #' @returns final_res a data frame of drug candidates prioritized by WCS, NCS, or Tau
  #' @author Kewalin Samart

  require(signatureSearch)
  # get CMap 2.0 connectivity scores
  qsig <- qSig(query = list(upset = up_genes,
                            downset = dn_genes),
               gess_method = "LINCS",
               refdb = db_path)
  # obtain query results
  res <- gess_lincs(qsig, sortby = "Tau", tau = TRUE)

  final_res <- finalize_result(result(res))

  return(final_res)
}

compute_Cor_based_scores <- function(full_signature_matrix, score_method, db_path){
  #' @description This function quantifies drug candidates for a given disgenes_values prioritized by Cor methods: XCor, XSpe
  #' @param full_signature_matrix a dataframe of expression values with associated disease genes as row names
  #' @param score_method types of correlatin e.g. "Cor_spearman", "Cor_pearson"
  #' @param db_path obtained form either setup_cmap1db() or setup_lincsdb() functions
  #' @returns final_res a data frame of drug candidates prioritized by XCor, XSpe
  #' @author Kewalin Samart

  require(signatureSearch)
  # get CMap 2.0 connectivity scores
  qsig <- qSig(query = full_signature_matrix,
               gess_method = "Cor",
               refdb = db_path)
  # obtain query results
  if(score_method == "Cor_spearman"){
    res <- gess_cor(qsig, method = "spearman")
  }else if(score_method == "Cor_pearson") {
    res <- gess_cor(qsig, method = "pearson")
  }

  final_res <- finalize_result(result(res))

  return(final_res)
}

combine_res_meta <- function(final_res){
  #' @description This function merge final_res with drug metadata
  #' @param final_res a data frame of drug candidates
  #' @returns final_res a data frame of drug candidates with drug info
  #' @author Kewalin Samart

  # get perturbation metdata
  pert_meta <- read.csv(file=here("data/metadata/repurposing_drugs_20200324.csv"),skip=9)
  pert_meta <- pert_meta[c("pert_iname","clinical_phase")]
  # merge metadata with the final result
  final_res <- merge(final_res, pert_meta, by.x = "pert", by.y = "pert_iname", all.y = FALSE)

  return(final_res)
}

get_FDAapproved_drugs <- function(final_res){
  #' @description This function filters out non FDA-approved drugs from final_res
  #' @param final_res a data frame of drug candidates with drug info
  #' @returns final_res a data frame of FDA-approved drug candidates with drug info
  #' @author Kewalin Samart

  # filter only FDA-approved drugs
  final_res <- final_res[final_res$clinical_phase == "Launched", ]

  # reorder the result based on the score_method
  # final_res <- final_res[order(abs(final_res[,score_method]), decreasing = TRUE),]
  # reset row names
  rownames(final_res) <- 1:nrow(final_res)

  return(final_res)
}

finalize_result <- function(final_res){
  #' @description This function finalize a data frame of drug candidates:
  #' @detail combine metadata, filter FDA-approved drugs, remove rows with NAs
  #' @param final_res a data frame of drug candidates
  #' @returns final_res a complete data frame of drug candidates
  #' @author Kewalin Samart

  # get perturbation metdata and merge with the final result
  final_res <- combine_res_meta(final_res)
  final_res <- get_FDAapproved_drugs(final_res)
  final_res <- na.omit(final_res)
  rownames(final_res) <- 1:nrow(final_res)

  return(final_res)
}
