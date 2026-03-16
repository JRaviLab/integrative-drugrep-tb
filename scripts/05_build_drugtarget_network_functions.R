# functions to build disease-drug target network given drug names, drug perturbed genes (optional), disease genes (optional), disease enriched pathways (optional)
# node: drugs, enriched pathways from Drug Set Enrichment Analysis (DSEA), drug-target genes, MOAs
# edges:  1. get drug-target information from  DrugBank, PharmGKB, ChEMBL, Drug Target Commons,and others (https://www.dgidb.org/browse/sources)
#         2. get drug-drug_gene information from LINCS level 5; significantly perturbed drug genes based on z-scores (both up- and downregulated)
#         3. get drug-moa information from CLUE via signatureSearch
#         4. get (drug) target-target connections from STRING to connect each MOA modules
#         5. get drug-pathway connections from DSEA via signatureSearch
# P.S. MOA-genes not included as drugrep hub data does not provide unique sets of MOA-genes for each MOA
#      I suppose they only give intersecting MOA-target-genes when a drug is associated with >1 MOAs (code at the bottom)

source(here("scripts/03_summarize_drugs_methodswise_functions.R"))
library(signatureSearch)

get_DTG_interactions_DGIdb <- function(drug_names,
                                      interaction_file = here("data/metadata/interactions_DGIdb.tsv"),
                                      values = "all",
                                      tail = "right",
                                      percentile = 0.9){
  #' @description this function obtains drug-target interactions from interactions file provided by DGIdb (https://www.dgidb.org/downloads)
  #' @param interaction_file path to interactions file from DGIdb
  #' @param values a string indicating which non-zero values to include in the distribution: "neg", "pos", if other values specified for this argument, it means including all values.
  #' @param tail a string indicating which side of the distribution to get the indicated top percentile from: "right", "left"
  #' @param percentile a numeric indicating top percentile scale for grabbing top drug-gene interactions; set to 0.9 by default
  #' @returns  DGIdb_targets_list a list of character vectors of drug-target interactions after filtering out interactions with scores lower than the given percentile threshold

  require(readr)

  print(" ======= Getting gene target information from DGIdb =======")
  if(tail == "left"){
    percentile = 1 - percentile
  }
  print(paste0("Threshold for drug-gene interactions  at ", percentile*100, "th percentile (getting top ",(percentile)*100,"% of the scores from the ",tail," tail)" ))

  DGIdb_interactions_df <- read_tsv(interaction_file)
  DGIdb_interactions_RA <- DGIdb_interactions_df[tolower(DGIdb_interactions_df$drug_name) %in% drug_names,]
  DGIdb_interactions_RA$drug_name <- tolower(DGIdb_interactions_RA$drug_name)

  DGIdb_targets_list <- list()
  for(drug in DGIdb_interactions_RA$drug_name){
    DGIdb_intRA_drug <- DGIdb_interactions_RA[DGIdb_interactions_RA$drug_name == drug,]
    drug_int_scores <- as.numeric(DGIdb_intRA_drug$interaction_score)
    drug_int_scores <- drug_int_scores[!is.na(drug_int_scores)]
    # grab the interactions within the top percentile
    intscore_threshold = determine_threshold(score_vec = drug_int_scores , score_percentile = percentile, values, tail)
    DGIdb_intRA_drug_toppct <- DGIdb_intRA_drug[DGIdb_intRA_drug$interaction_score > intscore_threshold,]
    DGIdb_top_targets <- unique(DGIdb_intRA_drug_toppct$gene_name)
    DGIdb_top_targets <- DGIdb_top_targets[!(DGIdb_top_targets == "NULL")]
    DGIdb_top_targets <- DGIdb_top_targets[!is.na(DGIdb_top_targets)]
    # compute averaged interaction score for each drug-target pair
    avg_intscore_vec <- c()
    for(target in DGIdb_top_targets){
      target_intscores <- DGIdb_intRA_drug_toppct[DGIdb_intRA_drug_toppct$gene_name == target,]$interaction_score
      target_intscores <- target_intscores[!is.na(target_intscores)]
      target_intscores <- target_intscores[!(target_intscores == "NULL")]
      target_intscores <- as.numeric(target_intscores)
      avg_target_intscore <- mean(target_intscores)
      avg_intscore_vec <- c(avg_intscore_vec,avg_target_intscore )
    }

    if(length(DGIdb_top_targets) != 0){
      names(avg_intscore_vec) <- DGIdb_top_targets
      DGIdb_targets_list[[drug]] <- avg_intscore_vec
    }
  }
  return(DGIdb_targets_list)
}

get_DGIdb_targets_table <- function(DGIdb_targets_list){
  #' @description this function convert list of drug targets (DGIdb_targets_list) to data frame
  #' @param DGIdb_targets_list a list of character vectors of drug-target interactions after filtering out interactions with scores lower than the given percentile threshold
  #' @returns DGIdb_targets_df a data frame with drug-target information; columns: "drug.target","drug_name","target","avg_interaction_score"

  require(tidyverse)
  DGIdb_targets_df <- as.data.frame(unlist(DGIdb_targets_list))
  DGIdb_targets_df$drug.target <- row.names(DGIdb_targets_df)
  colnames(DGIdb_targets_df)[1] <- "avg_interaction_score"
  DGIdb_targets_df <- DGIdb_targets_df %>% separate(remove=FALSE,col = drug.target,into=c("drug_name","target"),sep = "\\.")
  row.names(DGIdb_targets_df) <- NULL
  DGIdb_targets_df <- DGIdb_targets_df[c("drug.target","drug_name","target","avg_interaction_score")]

  return(DGIdb_targets_df)
}

get_drug_perturbed_genes <- function(drug_names,
                                     signature_matrix,
                                     signature_metadata,
                                     time_point = "24 h",
                                     dosage = "10 µM",
                                     N = 5,
                                     gene_conversion_file = here("data/metadata/Homo_sapiens.gene_info.tsv")){
  #' @description this function identifies significantly perturbed genes of given drug names based on LINCS level 5 signatures GSE70138
  #' @param drug_names a character vector containing drug names to identify enriched/perturbed drug genes
  #' @param signature_matrix a matrix of drug signatures where first col: rid represents gene ids and other columns are signature ids
  #' @example  signature_matrix = read_tsv("./data/drug_data/GSE92742_level5_signatures_lm.tsv")
  #' @param signature_metadata metadata table descriping details of the signature ids in the corresponding signature matrix
  #' @example signature_metadata = read_tsv("./data/metadata/GSE92742_Broad_LINCS_sig_info.txt")
  #' @param time_point a string indicating drug treatment duration
  #' @param dosage a string indicating drug concentration
  #' @param N a positive integer indicating number of top drug perturbed genes for both up- and downregulated
  #' @returns drug_gene_df a dataframe containing columns: "drug_name", "drug_gene"

  i = 1
  for(drug_name in drug_names){
    # obtain signatures ids of the given drug name
    signature_metadata = signature_metadata[which((signature_metadata$pert_itime == time_point) & (signature_metadata$pert_idose == dosage)),]
    drug_sigids = signature_metadata[which(signature_metadata$pert_iname == drug_name),]$sig_id
    # subset signature matrix by the drug signature ids
    row.names(signature_matrix) <- signature_matrix$rid
    drug_signature_matrix = signature_matrix[,which(colnames(signature_matrix) %in% drug_sigids)]
    row.names(drug_signature_matrix) <- signature_matrix$rid
    # compute average across cell lines
    drug_signature_matrix = as.matrix(drug_signature_matrix)
    avg_drug_sigmat = rowMeans(drug_signature_matrix)
    # sort by z-scores in descending order
    sorted_avg_drug_sigmat = avg_drug_sigmat[order(unname(avg_drug_sigmat), decreasing = TRUE)]
    # grab 20 genes from the head and tail of the drug-gene vector
    drug_perturbed_genes = names(c(head((sorted_avg_drug_sigmat),N),tail(sorted_avg_drug_sigmat,N)))
    drug_genes_df_ = as.data.frame(drug_perturbed_genes)
    colnames(drug_genes_df_)[1] <- "drug_gene"
    drug_genes_df_$drug_name <- drug_name

    if(i == 1){
      drug_genes_df = drug_genes_df_
    }else{
      drug_genes_df = rbind(drug_genes_df, drug_genes_df_)
    }
    i = i + 1
  }
  # convert GeneID to Symbol
  Homo_sapiens_gene_info <- read_tsv(gene_conversion_file)
  dp_genes_table <- Homo_sapiens_gene_info[which(Homo_sapiens_gene_info$GeneID %in% drug_genes_df$drug_gene),]
  drug_genes_symb_df <- merge(drug_genes_df, dp_genes_table, by.x = "drug_gene", by.y = "GeneID", all.x = TRUE, all.y = FALSE)
  drug_genes_df <- drug_genes_symb_df[c("drug_name","Symbol")]
  colnames(drug_genes_df)[2] <- "drug_gene"

  return(drug_genes_df)
}

get_drugtarget_NetworkNodeData <- function(drug_names,
                                          drug_genes_df,
                                          disease_pathways = c(),
                                          N = 15,
                                          interaction_file = here("data/metadata/interactions_DGIdb.tsv"),
                                          intscore_values = "all",
                                          intscore_tail = "right",
                                          intscore_pct = 0.9,
                                          moas = c()){
  #' @description this function takes drug names and build drug-target networks including pathway, MOAs,
  #' STRING gene-gene interaction information
  #' @param drug_names a character vector of drug names
  #' @param drug_genes_df a dataframe containing columns: "drug_name", "drug_gene"
  #' @param disease_pathways a character vector of disease pathways
  #' @param interaction_file path to interactions file from DGIdb
  #' @param intscore_values a string indicating which non-zero values to include in the distribution: "neg" (only negatives), "pos" (only positives), if other values specified for this argument, it means including all values.
  #' @param intscore_tail a string indicating which side of the distribution to get the indicated top percentile from: "right", "left"
  #' @param intscore_pct a numeric indicating top percentile scale for grabbing top drug-gene interactions; set to 0.9 by default
  #' @param moas a subset of moa to generate subnetwork of interest
  #' @returns drug_dpgene_tg_pw_moas a data frame containing node data for building a drug-target network with columns: "drug_name","moa","disease_pathway","enriched_pathway","target" or without "disease_pathway"

  require(signatureSearch)
  require(tidyverse)

  # get unique drug names
  drug_names <- unique(drug_names)

  # get list of drug target genes/proteins
  DGIdb_targets_list <- get_DTG_interactions_DGIdb(drug_names,
                                                   interaction_file,
                                                   values = intscore_values,
                                                   tail = intscore_tail,
                                                   percentile = intscore_pct)
  # create drug-target data frame
  DGIdb_targets_df <- get_DGIdb_targets_table(DGIdb_targets_list)

  print(" ======= Getting drug-enriched pathways using DSEA =======")
  # get drug-enriched pathways via DSEA (Drug Set Enrichment Analysis)
  hyperG_res <- dsea_hyperG(drugs=drug_names, type="GO", ont="BP")
  hyperG_res_df <- result(hyperG_res)
  # create drug-enriched pathways data frame
  drug_pathways_list <- list()
  for(drug in drug_names){
    drug_pathways_list_ <- list()
    for(i in 1:nrow(hyperG_res_df)){
      pathway <- hyperG_res_df$Description[i]
      GO_drugs_vec <- strsplit(hyperG_res_df$itemID[i],split = "/")
      if(drug %in% GO_drugs_vec[[1]]){
        drug_pathways_list_ <- append(drug_pathways_list_,pathway)
      }
      drug_pathways_list[[drug]] <- unique(drug_pathways_list_)
    }
  }
  drug_pathways_df <- as.data.frame(unlist(drug_pathways_list))
  colnames(drug_pathways_df)[1] <- "enriched_pathway"
  drug_pathways_df$drug_name <- gsub('[[:digit:]]+', '', row.names(drug_pathways_df))
  row.names(drug_pathways_df) <- NULL

  print(" ======= Getting drug-MOAs information from CLUE =======")
  # get MOAs information from CLUE for a given set of drug names
  data("clue_moa_list")
  drug_name_df <- as.data.frame(drug_names)
  colnames(drug_name_df)[1] <- "drug_name"
  drug_moas_df <- addMOA(as.data.frame(drug_name_df), "drug_name", clue_moa_list)
  # create drug-MOAs data frame
  drug_moas_list <- list()
  for(i in 1:nrow(drug_moas_df)){
    drug <- drug_moas_df$drug_name[i]
    moa_drugs_vec <- strsplit(drug_moas_df$MOAss[i],split = "; ")
    drug_moas_list[[drug]] <- moa_drugs_vec[[1]]
  }
  drug_moas_df <- as.data.frame(unlist(drug_moas_list))
  colnames(drug_moas_df)[1] <- "moa"
  drug_moas_df$drug_name <- gsub('[[:digit:]]+', '', row.names(drug_moas_df))
  row.names(drug_moas_df) <- NULL

  # build subnetwork based on the specified moas, if any
  if(length(moas) != 0){
    drug_moas_df = drug_moas_df[which(drug_moas_df$moa %in% moas),]
  }

  # get associated drug perturbed genes
  drug_dpgene <- merge(DGIdb_targets_df, drug_genes_df, by = "drug_name")
  # get associated moas for all drugs
  drug_dpgene_tg_pw <- merge(drug_dpgene, drug_pathways_df, by = "drug_name")
  # get associated known target genes for all drugs
  drug_dpgene_tg_pw_moas <-  merge(drug_moas_df, drug_dpgene_tg_pw, by = "drug_name")

  # add a column indicating whether the drug pathways are also one of the enriched disease pathways
  if(length(disease_pathways) != 0){
    drug_dpgene_tg_pw_moas = drug_dpgene_tg_pw_moas[which(drug_dpgene_tg_pw_moas$enriched_pathway %in% disease_pathways),]
  }

  drug_dpgene_tg_pw_moas <- drug_dpgene_tg_pw_moas[c("drug_name", "drug_gene","moa","enriched_pathway","target")]
  drug_dpgene_tg_pw_moas$drug.target <- NULL

  return(drug_dpgene_tg_pw_moas)
}

add_STRING_genes <- function(drug_dpgene_tg_pw_moas,
                             STRING_path = here("data/metadata/STRING_hs_interactions.txt"),
                             STRING_intscore = 0.7,
                             min_degree = 4,
                             disease_genes = c()) {
  #' @description this function adds PPIs connections from STRING to a subnetwork node table of interest drug_tg_pw_moas
  #' @param drug_dpgene_tg_pw_moas a data frame containing node data for building a drug-target network with columns: "drug_name", "drug_gene","moa","enriched_pathway","target"
  #' @param STRING_path path to human gene-gene interaction information file from STRING
  #' @param STRING_intscore a float indicating PPI interaction score cutoff to be included in the network; set to 0.7 by default
  #' @param disease_genes a character vector of disease-associated genes; these can be disease differentially expressed genes or other sources
  #' @returns network_nodes a data frame with columns: "drug_target1", "drug_name", "drug_gene", "moa", "enriched_pathway", "drug_target2", "tg_int_weight", "disease_target1",
  #' "disease_target2"
  #'
  print(" ======= Getting gene-gene interaction information from STRING =======")

  # load STRING interactions and filter by confidence score
  gene_gene_string <- read.delim(STRING_path, sep="\t", col.names = c("gene1","gene2","edge_weight"))
  gene_gene_string <- gene_gene_string[gene_gene_string$edge_weight >= STRING_intscore,]

  # identify network genes
  network_genes <- unique(union(drug_dpgene_tg_pw_moas$drug_gene, drug_dpgene_tg_pw_moas$target))

  # get only interactions where both genes are relevant
  gene_gene_string_subset <- gene_gene_string[
    (gene_gene_string$gene1 %in% network_genes | gene_gene_string$gene2 %in% network_genes), ]

  # identify additional STRING genes that are not in the original network
  all_genes <- unique(c(gene_gene_string_subset$gene1, gene_gene_string_subset$gene2))
  all_genes_butnw <- setdiff(all_genes, network_genes)

  # identify highly connected genes
  gene_gene_toadd <- data.frame()  # initialize empty dataframe
  for (gene in all_genes_butnw) {
    gene_gene_string_ <- gene_gene_string_subset[
      gene_gene_string_subset$gene1 %in% network_genes & gene_gene_string_subset$gene2 == gene, ]

    conn_num <- length(unique(gene_gene_string_$gene1))

    if (conn_num >= min_degree) {
      gene_gene_toadd <- rbind(gene_gene_toadd, gene_gene_string_)
    }
  }

  # merge new connections into the network
  network_nodes <- merge(drug_dpgene_tg_pw_moas, gene_gene_toadd,
                         by.x = "target", by.y = "gene1",
                         all.x = TRUE, all.y = FALSE)

  # rename columns
  colnames(network_nodes) <- gsub("target", "drug_target1", colnames(network_nodes))
  colnames(network_nodes) <- gsub("gene2", "drug_target2", colnames(network_nodes))
  colnames(network_nodes) <- gsub("edge_weight", "tg_int_weight", colnames(network_nodes))

  # identify disease genes
  network_nodes$disease_target1 <- as.integer(network_nodes$drug_target1 %in% disease_genes)
  network_nodes$disease_target2 <- as.integer(network_nodes$drug_target2 %in% disease_genes)

  # identify drug genes
  drug_genes <- unique(drug_dpgene_tg_pw_moas$drug_gene)
  network_nodes$drug_gene1 <- as.integer(network_nodes$drug_target1 %in% drug_genes)
  network_nodes$drug_gene2 <- as.integer(network_nodes$drug_target2 %in% drug_genes)

  return(network_nodes)
}

get_drugtarget_NetworkEdgeData <- function(network_nodes){
  #' @description this functions take node data and convert it to edge data with interaction types
  #' @param network_node_df a data frame consisting of the following column names:
  #' "drug_name","drug_gene","moa","enriched_pathway","drug_target1","drug_target2","edge_weight"
  #' @returns network_edges a data frame with "from", "to", and "type" indicating
  #' undirected and unweighted interactions with their interaction types for a given network node data frame.

  # target-target
  target_target <- network_nodes[c("drug_target1","drug_target2")]
  target_target <- unique(na.omit(target_target))
  colnames(target_target) <- c("from","to")
  target_target$type <- "target-target"

  # drug-dpgene
  drug_dpgene <- network_nodes[c("drug_name","drug_gene")]
  drug_dpgene <- unique(na.omit(drug_dpgene))
  colnames(drug_dpgene) <- c("from","to")
  drug_dpgene$type <- "drug-dpgene"

  # drug-target
  drug_target <- network_nodes[c("drug_name","drug_target1")]
  drug_target <- unique(na.omit(drug_target))
  colnames(drug_target) <- c("from","to")
  drug_target$type <- "drug-target"

  # drug-pathway
  drug_pathway <- network_nodes[c("drug_name","enriched_pathway")]
  drug_pathway <- unique(na.omit(drug_pathway))
  colnames(drug_pathway) <- c("from","to")
  drug_pathway$type <- "drug-pathway"

  # moa-drug
  drug_moa <- network_nodes[c("drug_name","moa")]
  drug_moa <- unique(na.omit(drug_moa))
  colnames(drug_moa) <- c("from","to")
  drug_moa$type <- "drug-moa"

  # combine drug-dpgene and drug-target and remove duplicates by flavoring drug-dpgene
  drug_gene <- unique(rbind(drug_dpgene, drug_target))
  drug_gene$dup_check <- paste0(drug_gene$from,"-",drug_gene$to)
  drug_gene <- drug_gene[!duplicated(drug_gene$dup_check), ]
  drug_gene$dup_check <- NULL

  # combine all edges
  network_edges <- unique(rbind(target_target, drug_gene, drug_pathway, drug_moa))

  return(network_edges)
}
