# functions to obtain background genes for Pathway analyses and suchs -- Homo Sapiens
# created date: 07/21/22
# last modified: 04/11/23
# Kewalin Samart

GO_genes <- function(){
  #' @description this function returns all human gene ids present in the GO db
  #' @param None
  #' @returns all GO human gene ids
  #' @source https://support.bioconductor.org/p/124061/
  #' @author Kewalin Samart
  require(org.Hs.eg.db)
  df <- as.data.frame(org.Hs.egGO)
  go_genes <-  unique(sort(df$gene_id))

  return(go_genes)
}

KEGG_genes <- function(){
  #' @description this function returns all human gene ids present in the KEGG db
  #' @param None
  #' @returns all KEGG human gene ids
  #' @source https://support.bioconductor.org/p/124061/
  #' @author Kewalin Samart
  require(org.Hs.eg.db)
  df = as.data.frame(org.Hs.egPATH)
  kegg_genes = unique(sort(df$gene_id))

  return(kegg_genes)
}

bg_LINCS_genes <- function(GeneType="landmark"){
  #' @description select background genes from LINCS
  #' @param GeneType a string or a list of strings indicating one or more types of LINCS gene: (i) landmark (by default) (ii) inferred (iii) best inferred (iv) not inferred (v) reference
  #' @returns bg_lincs_genes a character vector of the selected LINCS GeneID
  #' @author Kewalin Samart

  lincs_genes <- read.delim("../data/metadata/LINCSGeneSpaceSub.txt", sep="\t")

  if(GeneType=="reference"){ # all genes in LINCS
    selected_lincs_genes <- lincs_genes$`Entrez.ID`
  }else{
    selected_lincs_genes <- lincs_genes[lincs_genes$Type %in% c(GeneType),]$`Entrez.ID`
  }
  bg_lincs_genes <- as.character(selected_lincs_genes)

  return(bg_lincs_genes)
}

bg_from_data <- function(metadata_path, data_path, extension="") {
  #' @description Collect union of genes across all datasets (background genes)
  #' @param metadata_path Path to metadata file listing datasets
  #' @param data_path Path to DE result tables or signatures
  #' @param extension String like "up_none", "dn_landmark", "full_none" (optional)
  #' @return Character vector of unique genes across datasets
  #' @author Kewalin Samart

  all_genes <- character(0)

  data_to_run <- readr::read_tsv(metadata_path)
  data_to_run <- data_to_run[data_to_run$signature == 1,]

  for (i in seq_len(nrow(data_to_run))) {
    file_name <- data_to_run$SIGNATURE_NAME[i]
    file_path <- paste0(data_path,"/",file_name,extension)

    if (!file.exists(file_path)) next

    dataset <- readr::read_tsv(file_path) %>% tidyr::drop_na()

    all_genes <- c(all_genes, as.character(dataset$GeneID))
  }

  return(unique(all_genes))
}

get_bg_genes <- function(bg_source, metadata_path=NULL, data_path=NULL, extra_arg=NULL, extension=NULL){
  #' @description get background genes by the given data source name
  #' @param bg_source a string indicating the name of database: "GO", "KEGG", "LINCS"
  #' @param metadata_path path to metadata file describing DE results/signatures details
  #' @example "./data/metadata/"
  #' @param data_path path to DE data table/ signatures
  #' @example "./data/uniformly_processed/microarray/DEresults/", "./data/uniformly_processed/microarray/signatures/dn/"
  #' @param extra_arg a string or a list of strings indicating one or more types of LINCS gene: (i) landmark (by default) (ii) inferred (iii) best inferred (iv) not inferred (v) reference
  #' @param extension
  #' @returns bg_genes background genes (all genes present in the source database)
  #' @author Kewalin Samart
  if(bg_source == "LINCS"){
    if(is.null(extra_arg)){
      bg_genes <- bg_LINCS_genes()
    }else{
      bg_genes <- bg_LINCS_genes(GeneType=extra_arg)
    }
  }else if(bg_source == "KEGG"){
    bg_genes <- KEGG_genes()
  }else if(bg_source == "GO") {
    bg_genes <- GO_genes()
  }else if(bg_source == "input data") {
    bg_genes <- bg_from_data(metadata_path,data_path,extension = extension)
  }
  return(bg_genes)
}

final_bg_genes <- function(gene_set, bg_source, GeneType=NULL){
  #' @description Given a gene set of interest (GeneID i.e., Entrezid) and name of database to perform analyses, this function gives back the intersecting genes as the final background genes.
  #' @param gene_set a character vector of GeneIDs of interest
  #' @param bg_source a string indicating the name of database: "GO", "KEGG", "LINCS"
  #' @param GeneType (optional) a string indicating gene type(s) from the L1000 project: (i) "landmark" (by default) (ii) "inferred" (iii) "best inferred" (iv) "not inferred" (v) "reference"
  #' @returns final_bg_genes a character vector containing the final backgroud genes
  #' @author Kewalin Samart
  bg_genes = get_bg_genes(bg_source, extra_arg=GeneType)
  final_bg_genes = intersect(as.character(gene_set), bg_genes)

  return(final_bg_genes)
}