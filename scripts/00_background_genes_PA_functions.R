# functions to obtain background genes for Pathway analyses and suchs -- Homo Sapiens
# last modified: 06/26/25
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
  kegg_genes = as.character(unique(sort(df$gene_id)))

  return(kegg_genes)
}

bg_LINCS_genes <- function(GeneType="landmark"){
  #' @description select background genes from LINCS
  #' @param GeneType a string or a list of strings indicating one or more types of LINCS gene: (i) landmark (by default) (ii) inferred (iii) best inferred (iv) not inferred (v) reference
  #' @returns bg_lincs_genes a character vector of the selected LINCS GeneID
  #' @author Kewalin Samart
  require(here)
  lincs_genes <- read.delim(here("data/metadata/LINCSGeneSpaceSub.txt"), sep="\t")

  if(GeneType=="reference"){ # all genes in LINCS
    selected_lincs_genes <- lincs_genes$`Entrez.ID`
  }else{
    selected_lincs_genes <- lincs_genes[lincs_genes$Type %in% c(GeneType),]$`Entrez.ID`
  }
  bg_lincs_genes <- as.character(selected_lincs_genes)

  return(bg_lincs_genes)
}

bg_from_data <- function(metadata_path, data_path, extension=""){
  #' @description background genes across datasets
  #' @param metadata_path path to metadata file describing dataset details including SIGNATURE_NAME column
  #' with paths to signatures for background genes e.g., "data/v2/signatures/RNAseq_TB_signature_run_info.tsv"
  #' @param data_path path to DE data table/ signatures
  #' @example "data/v2/signatures/RNAseq/up",
  #' @param extension a string indicating a suffix of file names to get background genes from; "" by default; e.g., "up", "dn", "full"
  #' @returns genes all genes present across all input datasets
  #' @author Kewalin Samart

  union_bg_data_genes <- c()
  if(extension != ""){
    extension <- paste0("_",extension)
  }
  data_to_run <- read_tsv(here(metadata_path))
  for(i in 1:nrow(data_to_run)){
    # read in variables
    signature_name <- data_to_run$SIGNATURE_NAME[i]
    file_path <- here::here(data_path, paste0(signature_name, extension, ".tsv"))
      # check if file exists, if not then skip
      if(file.exists(file_path)){
        dataset <- read.delim(file_path, sep="\t")
      }else{
        next
      }
    }

    dataset <- na.omit(dataset) # remove rows with all NAs in the data

    if(i == 1){
      union_bg_data_genes = as.character(dataset$GeneID)
      print("yes, i is equal to 1")
    }else{
      union_bg_data_genes <- c(union_bg_data_genes,as.character(dataset$GeneID)) # get all the genes present across every input datasets
    }

  return(unique(union_bg_data_genes)) # return only unique genes
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
    bg_genes <- bg_from_data(metadata_path = metadata_path, data_path = data_path, extension = extension)
  }

  return(bg_genes)
}
