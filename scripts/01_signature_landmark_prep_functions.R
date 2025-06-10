# functions to filter genes with L1000 genes
# created date: 06/17/21
# latest modified: 04/09/24
# Kewalin Samart

# load the needed libraries
library(readr)
library(tidyverse)
library(dplyr)

getL1000 <- function(LINCSGenes_path = "./data/metadata/LINCSGeneSpaceSub.txt"){
  #' @description gets L1000 genes for filtering signatures
  #' @source: http://amp.pharm.mssm.edu/l1000fwd/download_page
  #' @param LINCSGenes_path path to file containing all the LINCS genes including landmarks, best inferred, and inferred
  #' @returns L1000_gene_symbols a list of L1000 gene symbols (978 in total)
  #' @author Kewalin Samart

  # read in L1000 gene file
  LINCSGenes <- read.delim(LINCSGenes_path, sep='\t')
  # get L1000 gene list
  L1000_genes <- LINCSGenes[,LINCSGenes$Type == 'landmark']
  L1000_gene_entrezid <- as.list(L1000_genes$`Entrez ID`)

  return(L1000_gene_entrezid)
}

prepare_signature <- function(signature, L1000 = FALSE){
  #' @description checks column names of the signature before proceeding down stream analyses
  #' @param signature
  #' @param L1000 a boolean indicating whether to suset only L1000 genes
  #' @returns signature with 'Ensembl', 'Symbol', 'GeneID', 'P.Value', 'adj.P.Val', 'log2FoldChange', and other columns of stats values
  #' @author Kewalin Samart

  # read in gene info for gene conversion
  gene_info <- read.delim(file = here("data/metadata/Homo_sapiens.gene_info.csv"))
  gene_info <- gene_info[,c("GeneID","Symbol","GeneID", "P.Value", "adj.P.Val")]
  # check if the input signature contains all the required column names
  required_colnames <- c('Ensembl', 'Symbol', 'GeneID', 'P.Value', 'adj.P.Val')
  sig_cols <- colnames(signature)

  if(required_colnames %in% sig_cols){
    print("The signature contains all the required column names: Ensembl, Symbol, GeneID, P.Value, adj.P.Val")
    return(signature)
  }else{
    missing_colnames <- required_colnames[which(!required_colnames %in% sig_cols)]
    print(paste0("The signature is missing ", missing_colnames))
  }
  # if missing Symbol...
  if("Symbol" %in%  missing_colnames){
    if("GeneID" %in% sig_cols){
      # convert GeneID (Entrezid) to Symbol
      signature <- merge.data.frame(signature,gene_info,by = "GeneID")
    }else if("Ensembl" %in% sig_cols){
      # convert Ensembl to Symbol
      signature <- merge.data.frame(signature,gene_info,by = "Ensembl")
    }
  }
  # if missing GeneID...
  if("GeneID" %in%  missing_colnames){
    if("Symbol" %in% sig_cols){
      # convert GeneID (Entrezid) to Symbol
      signature <- merge.data.frame(signature,gene_info,by = "Symbol")
    }else if("Ensembl" %in% sig_cols){
      # convert Ensembl to Symbol
      signature <- merge.data.frame(signature,gene_info,by = "Ensembl")
    }
  }
  # if missing Ensembl...
  if("Ensembl" %in%  missing_colnames){
    if("Symbol" %in% sig_cols){
      # convert GeneID (Entrezid) to Symbol
      signature <- merge.data.frame(signature,gene_info,by = "Symbol")
    }else if("GeneID" %in% sig_cols){
      # convert Ensembl to Symbol
      signature <- merge.data.frame(signature,gene_info,by = "GeneID")
    }
  }
  # finalize the processed signature (output)
  signature <- signature[,c("Ensembl", "GeneID","Symbol","log2FoldChange", "P.Value", "adj.P.Val")]
  processed_signature <- signature %>% arrange(adj.P.Val)

  if(L1000){
    processed_signature = filterSig_withL1000(processed_signature, L1000_genes=getL1000())
  }

  return(processed_signature)
}

filterSig_withL1000 <- function(signature, L1000_genes=getL1000()){
  #' @description filters signature with L1000 genes
  #' @param signature output from get_signature()
  #' @param L1000_genes list of L1000 gene ids; default val set to the output from getL1000()
  #' @returns L1000_filtered_sig signature with only L1000 genes
  #' @author Kewalin Samart

  # filter the signature dataframe with L1000 genes
  L1000_filtered_sig <- signature[which(signature$GeneID %in% L1000_genes),]

  return(L1000_filtered_sig)
}
