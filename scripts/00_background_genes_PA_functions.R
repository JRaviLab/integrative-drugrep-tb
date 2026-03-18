# functions to obtain background genes for Pathway analyses -- Homo Sapiens

library(tidyverse)

GO_genes <- function() {
  #' @description this function returns all human gene ids present in the GO db
  #' @param None
  #' @returns all GO human gene ids
  #' @source https://support.bioconductor.org/p/124061/
  require(org.Hs.eg.db)
  df <- as.data.frame(org.Hs.egGO)
  go_genes <- unique(sort(df$gene_id))

  return(go_genes)
}

KEGG_genes <- function() {
  #' @description this function returns all human gene ids present in the KEGG db
  #' @param None
  #' @returns all KEGG human gene ids
  #' @source https://support.bioconductor.org/p/124061/
  require(org.Hs.eg.db)
  df <- as.data.frame(org.Hs.egPATH)
  kegg_genes <- unique(sort(df$gene_id))

  return(kegg_genes)
}

bg_LINCS_genes <- function(GeneType = "landmark") {
  #' @description select background genes from LINCS
  #' @param GeneType a string or a list of strings indicating one or more types of LINCS gene: (i) landmark (by default) (ii) inferred (iii) best inferred (iv) not inferred (v) reference
  #' @returns bg_lincs_genes a character vector of the selected LINCS GeneID

  lincs_genes <- read.delim(here("data/metadata/LINCSGeneSpaceSub.txt"), sep = "\t")

  if (GeneType == "reference") { # all genes in LINCS
    selected_lincs_genes <- lincs_genes$`Entrez.ID`
  } else {
    selected_lincs_genes <- lincs_genes[lincs_genes$Type %in% c(GeneType), ]$`Entrez.ID`
  }
  bg_lincs_genes <- as.character(selected_lincs_genes)

  return(bg_lincs_genes)
}

bg_from_data <- function(metadata_path, data_path, direction, extension = ".tsv") {
  #' @description Collect union of genes across all datasets (background genes)
  #' @param metadata_path Path to metadata file listing datasets
  #' @param data_path Path to DE result tables or signatures
  #' @param direction options: "up", "dn"
  #' @param extension String like "_up.tsv", "_dn.tsv" (optional)
  #' @return Character vector of unique genes across datasets

  all_genes <- character(0)

  data_to_run <- readr::read_tsv(metadata_path)
  data_to_run <- data_to_run[data_to_run$signature == 1, ]

  for (i in seq_len(nrow(data_to_run))) {
    if(grepl("DE_results", data_path)){
      file_name <- data_to_run$SIGNATURE_FULL_NAME[i]
    }else{
      file_name <- data_to_run$SIGNATURE_NAME[i]
    }
    file_path <- paste0(data_path, "/", direction, "/", file_name, extension)
    print(file_path)

    if (!file.exists(file_path)) next

    dataset <- readr::read_tsv(file_path) %>% tidyr::drop_na()

    all_genes <- c(all_genes, as.character(dataset$GeneID))
  }

  return(unique(all_genes))
}

get_bg_genes <- function(bg_source, metadata_path = NULL, data_path = NULL, lincs_genetype = NULL, extension = ".tsv", direction_input_data = NULL) {
  #' @description get background genes by the given data source name
  #' @param bg_source a string indicating the name of database: "GO", "KEGG", "LINCS"
  #' @param metadata_path path to metadata file describing DE results/signatures details
  #' @example "data/signatures/RNAseq_TB_signature_run_info.tsv"
  #' @param data_path path to DE data table/ signatures
  #' @example "data/DE_results/RNAseq", "data/DE_results/microarray"
  #' @param lincs_genetype a string or a list of strings indicating one or more types of LINCS gene: (i) landmark (by default) (ii) inferred (iii) best inferred (iv) not inferred (v) reference
  #' @param extension a string indicating file extension e.g., "_up.tsv", "_dn.tsv"
  #' @param direction_input_data options: "up", "dn"
  #' @returns bg_genes background genes (all genes present in the source database)
  if (bg_source == "LINCS") {
    if (is.null(lincs_genetype)) {
      bg_genes <- bg_LINCS_genes()
    } else {
      bg_genes <- bg_LINCS_genes(GeneType = lincs_genetype)
    }
  } else if (bg_source == "KEGG") {
    bg_genes <- KEGG_genes()
  } else if (bg_source == "GO") {
    bg_genes <- GO_genes()
  } else if (bg_source == "input data") {
    if (is.null(direction_input_data)) {
      stop("need to specify a valid direction if using input data as background genes")
    }
    bg_genes <- bg_from_data(metadata_path, data_path, direction = direction_input_data, extension = extension)
  }
  return(bg_genes)
}

final_bg_genes <- function(gene_set, bg_source, GeneType = NULL) {
  #' @description Given a gene set of interest (GeneID i.e., Entrezid) and name of database to perform analyses, this function gives back the intersecting genes as the final background genes.
  #' @param gene_set a character vector of GeneIDs of interest
  #' @param bg_source a string indicating the name of database: "GO", "KEGG", "LINCS"
  #' @param GeneType (optional) a string indicating gene type(s) from the L1000 project: (i) "landmark" (by default) (ii) "inferred" (iii) "best inferred" (iv) "not inferred" (v) "reference"
  #' @returns final_bg_genes a character vector containing the final backgroud genes
  bg_genes <- get_bg_genes(bg_source, lincs_genetype = GeneType)
  final_bg_genes <- intersect(as.character(gene_set), bg_genes)

  return(final_bg_genes)
}

get_combined_bg_genes <- function(metadata_path_rnaseq, data_path_rnaseq,
                                  metadata_path_microarray, data_path_microarray,
                                  direction, GO_genes_vector, combine) {
  #' @param metadata_path_rnaseq: metadata tsv for RNAseq
  #' @param data_path_rnaseq: path to DE results for RNAseq
  #' @param metadata_path_microarray: metadata tsv for microarray
  #' @param data_path_microarray: path to DE results for microarray
  #' @param direction: "up" or "dn"
  #' @param GO_genes_vector: vector of GO-annotated gene IDs
  #' @param combine: "intersect", "union"

  extension_rnaseq <- paste0("_", direction, ".tsv")
  extension_microarray <- paste0("_", direction, ".tsv")

  # RNAseq background
  bg_genes_rnaseq_raw <- get_bg_genes(
    bg_source = "input data",
    metadata_path = metadata_path_rnaseq,
    data_path = paste0(data_path_rnaseq),
    lincs_genetype = NULL,
    extension = extension_rnaseq,
    direction_input_data = direction
  )
  bg_genes_rnaseq <- intersect(bg_genes_rnaseq_raw, GO_genes_vector)

  # Microarray background
  bg_genes_microarray_raw <- get_bg_genes(
    bg_source = "input data",
    metadata_path = metadata_path_microarray,
    data_path = paste0(data_path_microarray),
    lincs_genetype = NULL,
    extension = extension_microarray,
    direction_input_data = direction
  )
  bg_genes_microarray <- intersect(bg_genes_microarray_raw, GO_genes_vector)

  if (combine == "intersect") {
    # Combine across technologies (intersection)
    bg_genes_combined <- intersect(bg_genes_rnaseq, bg_genes_microarray)
  } else if (combine == "union") {
    bg_genes_combined <- union(bg_genes_rnaseq, bg_genes_microarray)
  }

  return(bg_genes_combined)
}
