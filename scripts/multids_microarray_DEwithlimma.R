# script to perform differential gene expression on processed microarray data matrices
# using limma package: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
# ref (design matrix): https://rpubs.com/ge600/limma
# ref (microarray DE analysis on 2 groups): https://alexslemonade.github.io/refinebio-examples/02-microarray/differential-expression_microarray_01_2-groups.html
# last modified: 06/26/25
# author: Kewalin Samart

library(limma)
library(magrittr)
library(dplyr)
library(readr)
library(stringr)
library(here)

## set up arguments
args <- commandArgs(TRUE)
# args[1]: path to tsv file with datasets with arguments to perform DE analyses e.g., "./data/metadata/TB_microarray_args.tsv"
# args[2]: whether or not to filter DE results by "landmark", "best inferred", "inferred", "reference" (all gene types combined) genes in LINCS database, or "none"
# args[3]: adj.p val cutoff (0.05 by default)
args[1] <- "data/v2/microarray_data_forDE/clean_TB_sample_metadata_classification.tsv"
args[2] <- "landmark"
args[3] <- 0.05

# read in argument file
args_file = read.delim(here(args[1]), sep='\t')
lincs_genes = read.delim(here("data/metadata/LINCSGeneSpaceSub.txt"), sep='\t')
GeneType <- as.character(str_split(args[2],",")[[1]])

signature_boolean <- list()
platform_list <- list()

### modified up to here*****

for(i in 1:nrow(args_file)){
  meta_path = paste0(args_file$meta_path[i]) # e.g., "/data/metadata/metadata_GSE16250.tsv"
  exprmat_path = paste0(args_file$exprmat_path[i]) # e.g., "/data/expression_matrices/GSE16250.tsv"

  ## read in input files
  # read in metadata
  metadata <- read.delim(meta_path,sep = "\t")
  # read in expression matrix
  expr_mat <- read.delim(exprmat_path,sep = "\t")

  # get common sample ids
  common_sampleids = intersect(colnames(expr_mat),metadata$refinebio_accession_code)
  expression_df <- expr_mat[,common_sampleids]
  row.names(expression_df) <- expr_mat$Gene
  metadata <- metadata %>%filter(metadata$refinebio_accession_code %in% common_sampleids)

  ## pre-processing inputs
  # set up conditions for comparison/contrast
  condition_colname = args_file$condition_colname[i]
  inf_keywords = args_file$inf_keywords[i]
  control_keywords = args_file$control_keywords[i]

  inf_indices = grep(inf_keywords, metadata[,c(condition_colname)])
  inf_samples = metadata[c(inf_indices),]$refinebio_accession_code

  control_indices = grep(control_keywords, metadata[,c(condition_colname)])
  control_samples = metadata[c(control_indices),]$refinebio_accession_code

  # Check if this is in the same order # this should be done before specifying conditions
  expression_df = expression_df[metadata$refinebio_accession_code]
  all.equal(colnames(expression_df), metadata$refinebio_accession_code)

  # subset metadata and expression matrix by samples: infected and control samples
  selected_samples = c(inf_samples,control_samples)
  metadata <- metadata %>%filter(metadata$refinebio_accession_code %in% selected_samples)
  conditions = as.character(metadata[,c(condition_colname)]) # get condition names for DE comparison
  platform_ids = unique(metadata$platform_id) # get platform id(s) associated with the selected samples
  expression_df <-  expression_df[,selected_samples]

  # replace inf_indices with 'infected' and control_indices with 'control' # ** something is wrong here
  # the indices don't match
  # get new indices
  inf_indices = grep(inf_keywords, conditions)
  control_indices = grep(control_keywords, conditions)

  conditions[inf_indices ] <- "infected"
  conditions[control_indices ] <- "control"

  # create a design matrix
  groups <- factor(conditions, levels = c("infected","control"))

  if(length(platform_ids) > 1) {
    print("removing batch effect caused by platform")
    design_mat = model.matrix(~0 + groups + platform_ids)
    platform_ids <- cat(paste(platform_ids, collapse = "_"))
  }else if(length(platform_ids) == 1){
    design_mat = model.matrix(~0 + groups)
    print("no batch effect caused by GSE or platform")
  }

  # record platform ids
  platform_list[[i]] <- platform_ids

  colnames(design_mat) <- c("infected","control")

  ## perform differential gene expression analysis
  # apply linear model to the expression matrix
  fit <- lmFit(expression_df, design = design_mat)
  # compute batch-corrected mean log expression
  mean_exp = as.data.frame(fit$coefficients)

  # check and create mean_exp dirname
  mean_exp_dirname <- paste0("./batch_corrected_mean_exp/")
  if(!dir.exists(mean_exp_dirname)) {
    dir.create(mean_exp_dirname)
  }

  write_tsv(mean_exp, file = paste0(mean_exp_dirname, args_file$refine.bio_accession[i],"_",platform_ids,"_",args_file$file_name[i],"_batch_corrected_mean_log_expression.tsv"))

  # apply empirical Bayes to smooth standard errors
  fit <- eBayes(fit)

  # correct biases by performing multiple testing, and obtain result table
  stats_df <- topTable(fit, number = nrow(expression_df),) %>% tibble::rownames_to_column("Gene")

  ## add gene annotations
  # get gene annotation table
  annot_path = paste0("./annotation/Homo_sapiens.gene_info.csv")
  annotdata <- read.delim(annot_path,sep = ",")  # GeneID, Symbol, Ensembl
  ensembl <- stats_df$Gene
  annotdata_subset <- annotdata %>%filter(annotdata$Ensembl %in% ensembl)
  annotdata_subset = annotdata_subset[,c('GeneID','Symbol','Ensembl')]
  res_df =  merge(stats_df, annotdata_subset,  by.x = 'Gene', by.y = 'Ensembl', all.x = TRUE, all.y = FALSE)
  res_df = res_df[c("Gene", "GeneID","Symbol","infected","control","AveExpr","F","P.Value","adj.P.Val")]
  colnames(res_df)[1] <- "Ensembl"
  res_df$log2FoldChange = log2(res_df$infected) - log2(res_df$control)

  # identify duplicated genes (EntrezID)
  duplicates_booleans <- duplicated(res_df$GeneID)
  duplicated_genes <- res_df$GeneID[duplicates_booleans]
  duplicated_genes <- duplicated_genes[!is.na(duplicated_genes)]

  if(length(duplicated_genes) != 0){
    for(gene in duplicated_genes){
      # look for its maximum logFC value and keep the max
      gene_indices <- which(res_df$GeneID == gene)
      logFC_vec <- res_df$log2FoldChange[gene_indices]
      names(logFC_vec) <- gene_indices
      max_index <- as.numeric(names(logFC_vec)[which.max(abs(logFC_vec))])
      not_max_index <- gene_indices[which(gene_indices != max_index)]
      res_df <- res_df[-c(not_max_index), ]
    }
  }

  print(paste0("Saving DE results for Contrast: ",inf_keywords," vs ",control_keywords))

  # check and create DE dirname
  DE_dirname <- "./DE_results/"
  if(!dir.exists(DE_dirname)) {
    dir.create(DE_dirname)
  }

  ## save DE stats table
  write_tsv(res_df, file = paste0(DE_dirname,args_file$refine.bio_accession[i],"_",platform_ids,"_",args_file$file_name[i],"_",GeneType,".tsv"))

  if(!is.null(args[2])){
    if(args[2] == "none"){
    }else if(args[2] == "reference"){
      reference_genes <- as.character(lincs_genes$`Entrez.ID`)
      res_df <- res_df[as.character(res_df$GeneID) %in% reference_genes,]
    }else{
      subset_genes <- as.character(lincs_genes[lincs_genes$Type %in% GeneType,]$`Entrez.ID`)
      res_df <- res_df[as.character(res_df$GeneID) %in% subset_genes,]
    }

    # concatenate GeneType strings for file naming
    if(length(GeneType) > 1){
      GeneType <- paste(GeneType, collapse = "-")
    }else{
      GeneType <-toString(GeneType[1])
    }

    # filter out genes with adj.p val > 0.05 (not statistically significant)
    if(!is.null(args[3])){
      signi_cutoff <- as.numeric(args[3])
      res_df <- res_df[res_df$adj.P.Val < signi_cutoff,]
    }else{res_df <- res_df[res_df$adj.P.Val < 0.05,]}

    # check dimension of the signature based on the specified significance-level cutoff
    if(dim(res_df)[1] == 0){
      print(paste0('There is no significantly expressed genes at the specified adjusted p-value of ', signi_cutoff))
      print(paste(args_file$refine.bio_accession[i],"for",inf_keywords,"vs",control_keywords,sep=" "))
      signature_boolean[[i]] <- 0
      next
    }else{
      signature_boolean[[i]] <- 1
    }

    print(paste0("Saving signture for Contrast: ",inf_keywords," vs ",control_keywords))

    # check and create signature dirnames for up, dn, and full signatures
    sig_dirname <- "./signatures/"
    if(!dir.exists(sig_dirname)) {
      dir.create(sig_dirname)
    }
    sig_up_dirname <- "./signatures/up/"
    if(!dir.exists(sig_up_dirname)) {
      dir.create(sig_up_dirname)
    }
    sig_dn_dirname <- "./signatures/dn/"
    if(!dir.exists(sig_dn_dirname)) {
      dir.create(sig_dn_dirname)
    }
    sig_full_dirname <- "./signatures/full/"
    if(!dir.exists(sig_full_dirname)) {
      dir.create(sig_full_dirname)
    }

    ## save full signature
    write_tsv(res_df, file = paste0(sig_full_dirname,args_file$refine.bio_accession[i],"_",platform_ids,"_",args_file$file_name[i],"_full","_",GeneType,".tsv"))

    ## save up/dn signatures
    res_df_up <- res_df[res_df$log2FoldChange > 0,]
    res_df_dn <- res_df[res_df$log2FoldChange < 0,]

    print(paste0("Saving up/dn signtures for Contrast: ",inf_keywords," vs ",control_keywords))
    write_tsv(res_df_up, file = paste0(sig_up_dirname,args_file$refine.bio_accession[i],"_",platform_ids,"_",args_file$file_name[i],"_up","_",GeneType,".tsv"))
    write_tsv(res_df_dn, file = paste0(sig_dn_dirname,args_file$refine.bio_accession[i],"_",platform_ids,"_",args_file$file_name[i],"_dn","_",GeneType,".tsv"))
  }
}

# saving signature generation record
args_file$platform <- as.character(unlist(platform_list))
args_file$signature <- as.numeric(unlist(signature_boolean))
date_time <- gsub(" ", "_",  Sys.time())
write_tsv(args_file, file = paste0(sig_dirname,date_time,"_marray_signatures_info.tsv"))
