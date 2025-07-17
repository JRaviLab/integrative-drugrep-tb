# script to perform pathway analyses: GO/KEGG/Reactome ORA
# last modified: 07/16/25
# Kewalin Samart

library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(readr)
library(here)

source(here("scripts/00_background_genes_PA_functions.R"))

## set up variables

technologies = c("microarray","RNAseq")

metadata_path = here(paste0("data/v2/signatures/",technology,"_TB_signature_run_info.tsv"))
signature_info = read_tsv(metadata_path)
signature_info = signature_info[signature_info$signature == 1,]

# path to DE data table/ signatures
data_path <- here(paste0("data/v2/DE_results/", technology))

# get file names
filenames <- list.files(path = data_path, all.files = TRUE, pattern = "*.tsv",
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
signature_filenames <- filenames[grepl(paste(signature_info$series_id, collapse = "|"), filenames)]

##### get up/dn DE results #####
library(tidyverse)
library(here)

technologies <- c("microarray", "RNAseq")

for (technology in technologies) {

  # Load signature metadata
  metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
  signature_info <- read_tsv(metadata_path) %>%
    filter(signature == 1)

  # Paths
  data_path <- here(paste0("data/v2/DE_results/", technology))
  up_folder <- file.path(data_path, "up")
  dn_folder <- file.path(data_path, "dn")

  # Create up/ and dn/ folders if missing
  dir.create(up_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(dn_folder, showWarnings = FALSE, recursive = TRUE)

  # Get filenames matching series_ids
  filenames <- list.files(path = data_path, pattern = "\\.tsv$", full.names = FALSE)
  signature_filenames <- filenames[grepl(paste(signature_info$series_id, collapse = "|"), filenames)]

  # Process each signature file
  for (file in signature_filenames) {

    # Read DE table
    file_path <- file.path(data_path, file)
    de_df <- read_tsv(file_path)

    # Filter upregulated and downregulated
    up_genes <- de_df %>% filter(log2FoldChange > 0)
    dn_genes <- de_df %>% filter(log2FoldChange < 0)

    # Create new filenames with suffixes
    base_filename <- sub("\\.tsv$", "", file)
    up_file <- file.path(up_folder, paste0(base_filename, "_up.tsv"))
    dn_file <- file.path(dn_folder, paste0(base_filename, "_dn.tsv"))

    # Save
    write_tsv(up_genes, up_file)
    write_tsv(dn_genes, dn_file)
  }
}

##########

##### Pathway enrichment analysis #####

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(here)

#--------- Set up background genes ---------

get_combined_bg_genes <- function(metadata_path_rnaseq, data_path_rnaseq,
                                  metadata_path_microarray, data_path_microarray,
                                  direction, GO_genes_vector, combine) {
  #' @param metadata_path_rnaseq: Metadata TSV for RNAseq
  #' @param data_path_rnaseq: Path to DE results for RNAseq
  #' @param metadata_path_microarray: Metadata TSV for microarray
  #' @param data_path_microarray: Path to DE results for microarray
  #' @param direction: "up" or "dn"
  #' @param GO_genes_vector: Vector of GO-annotated gene IDs

  # Build extension
  extension_rnaseq <- paste0("_DESeq2_", direction, ".tsv")
  extension_microarray <- paste0("_limma_", direction, ".tsv")

  # RNAseq background
  bg_genes_rnaseq_raw <- get_bg_genes(
    bg_source = "input data",
    metadata_path = metadata_path_rnaseq,
    data_path = paste0(data_path_rnaseq,"/",direction),
    extension = extension_rnaseq
  )
  bg_genes_rnaseq <- intersect(bg_genes_rnaseq_raw, GO_genes_vector)

  # Microarray background
  bg_genes_microarray_raw <- get_bg_genes(
    bg_source = "input data",
    metadata_path = metadata_path_microarray,
    data_path = paste0(data_path_microarray,"/",direction),
    extension = extension_microarray
  )
  bg_genes_microarray <- intersect(bg_genes_microarray_raw, GO_genes_vector)

  if(combine == "intersect"){
    # Combine across technologies (intersection)
    bg_genes_combined <- intersect(bg_genes_rnaseq, bg_genes_microarray)
  }else if(combine == "union"){
    bg_genes_combined <- union(bg_genes_rnaseq, bg_genes_microarray)
  }


  return(bg_genes_combined)
}

# Define input paths
metadata_path_rnaseq <- here("data/v2/signatures/RNAseq_TB_signature_run_info.tsv")
data_path_rnaseq <- here("data/v2/DE_results/RNAseq")

metadata_path_microarray <- here("data/v2/signatures/microarray_TB_signature_run_info.tsv")
data_path_microarray <- here("data/v2/DE_results/microarray")

# Load GO genes
GO_genes_vector <- GO_genes()

# Get upregulated combined background
bg_genes_combined_up <- get_combined_bg_genes(
  metadata_path_rnaseq, data_path_rnaseq,
  metadata_path_microarray, data_path_microarray,
  direction = "up",
  GO_genes_vector,
  combine = "union"
)

# Get downregulated combined background
bg_genes_combined_dn <- get_combined_bg_genes(
  metadata_path_rnaseq, data_path_rnaseq,
  metadata_path_microarray, data_path_microarray,
  direction = "dn",
  GO_genes_vector,
  combine = "union"
)


#-------------------------------------------

technology <- "RNAseq"
directions <- c("up", "dn")
metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
data_path <- here(paste0("data/v2/DE_results/", technology))

for (direction in directions) {
  print(direction)
  # Set output directory
  dirname <- here(paste0("data/v2/", technology, "_pathways/", direction))
  if (!dir.exists(dirname)) {
    dir.create(dirname, recursive = TRUE)
  }
  dirname_read <- file.path(data_path, direction)
  filenames <- list.files(path = dirname_read, pattern = "*.tsv", full.names = FALSE)

  # Iterate through each file
  itr <- 1
  for (file in filenames) {

    print(paste0("Running iteration: ", itr))
    itr <- itr + 1

    genes_df <- read_tsv(file.path(dirname_read, file))

    # Use gene IDs as input for ORA
    gene_vector <- genes_df$log2FoldChange
    names(gene_vector) <- as.character(genes_df$GeneID)

    gene_list <- names(gene_vector)

    # Run ORA
    print(paste0("Getting enriched GO pathways for ", file))
    if(direction == "up"){
      bg_genes = bg_genes_combined_up
    }else if(direction == "dn"){
      bg_genes = bg_genes_combined_dn
    }

    enrichGO_res <- enrichGO(gene = gene_list,
                             OrgDb = org.Hs.eg.db,
                             readable = TRUE,
                             ont = "BP",
                             universe = bg_genes,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.1,
                             minGSSize = 5,
                             maxGSSize = 200)

    enrichGO_res <- as.data.frame(enrichGO_res)

    # Save results
    base_filename <- sub(paste0("\\.", direction, "\\.tsv$"), "", file)
    output_file <- file.path(dirname, paste0("GO_ORA_BGcorrected_", base_filename))

    write_tsv(enrichGO_res, file = output_file)
  }
}

###--------------- Aggregated signatures
bg_genes_source <- "input data"

### RNAseq
technology <- "RNAseq"
metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
data_path <- here(paste0("data/v2/DE_results/", technology))
if(technology == "RNAseq"){extension <- paste0("_DESeq2_",direction,".tsv")}
if(technology == "microarray"){extension <- paste0("_limma_",direction,".tsv")}
### up signatures
direction <- "up"
# Define background genes
dirname_read <- file.path(data_path, direction)
#up_aggr_signature_rnaseq <- as.character(read_tsv(here("data/v2/signatures/RNAseq/aggregated_signatures/up_aggregated_signature.tsv"))$GeneID)
up_aggr_signature_rnaseq <- as.character(read_tsv(here("data/v2/signatures/RNAseq/aggregated_signatures/RNAseq_up_aggregated_signature_0.8.tsv"))$GeneID)
enrichGO_res_up_aggr_rnaseq <- enrichGO(gene = up_aggr_signature_rnaseq,
                         OrgDb = org.Hs.eg.db,
                         readable = TRUE,
                         ont = "BP",
                         universe = bg_genes_combined_up,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         minGSSize = 5,
                         maxGSSize = 200)

enrichGO_res_up_aggr_rnaseq <- as.data.frame(enrichGO_res_up_aggr_rnaseq)

# Save results
output_file <- file.path(here("data/v2/RNAseq_pathways/up"), paste0("GO_ORA_BGcorrected_up_aggregated_RNAseq_0.8.tsv"))

write_tsv(enrichGO_res_up_aggr_rnaseq, file = output_file)


### dn signatures
direction <- "dn"
if(technology == "RNAseq"){extension <- paste0("_DESeq2_",direction,".tsv")}
if(technology == "microarray"){extension <- paste0("_limma_",direction,".tsv")}
# Define background genes
dirname_read <- file.path(data_path, direction)
dn_aggr_signature_rnaseq <- as.character(read_tsv(here("data/v2/signatures/RNAseq/aggregated_signatures/RNAseq_dn_aggregated_signature_0.8.tsv"))$GeneID)
enrichGO_res_dn_aggr_rnaseq <- enrichGO(gene = dn_aggr_signature_rnaseq,
                                        OrgDb = org.Hs.eg.db,
                                        readable = TRUE,
                                        ont = "BP",
                                        universe = bg_genes_combined_dn,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.1,
                                        minGSSize = 5,
                                        maxGSSize = 200)

enrichGO_res_dn_aggr_rnaseq <- as.data.frame(enrichGO_res_dn_aggr_rnaseq)

# Save results
output_file_dn_aggr_rnaseq <- file.path(here("data/v2/RNAseq_pathways/dn"), paste0("GO_ORA_BGcorrected_dn_aggregated_RNAseq_0.8.tsv"))

write_tsv(enrichGO_res_up_aggr_rnaseq, file = output_file_dn_aggr_rnaseq )

### microarray
technology <- "microarray"
metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
data_path <- here(paste0("data/v2/DE_results/", technology))
### up signatures
direction <- "up"
if(technology == "RNAseq"){extension <- paste0("_DESeq2_",direction,".tsv")}
if(technology == "microarray"){extension <- paste0("_limma_",direction,".tsv")}
# Define background genes
dirname_read <- file.path(data_path, direction)
up_aggr_signature_marray <- as.character(read_tsv(here("data/v2/signatures/microarray/aggregated_signatures/microarray_up_aggregated_signature_0.8.tsv"))$GeneID)
enrichGO_res_up_aggr_marray <- enrichGO(gene = up_aggr_signature_marray,
                                        OrgDb = org.Hs.eg.db,
                                        readable = TRUE,
                                        ont = "BP",
                                        universe = bg_genes_combined_up,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.1,
                                        minGSSize = 5,
                                        maxGSSize = 200)

enrichGO_res_up_aggr_marray <- as.data.frame(enrichGO_res_up_aggr_marray)

# Save results
output_file_up_aggr_marray <- file.path(here("data/v2/microarray_pathways/up"), paste0("GO_ORA_BGcorrected_up_aggregated_microarray_0.8.tsv"))

write_tsv(enrichGO_res_up_aggr_marray, file = output_file_up_aggr_marray)

### dn signatures
direction <- "dn"
if(technology == "RNAseq"){extension <- paste0("_DESeq2_",direction,".tsv")}
if(technology == "microarray"){extension <- paste0("_limma_",direction,".tsv")}
# Define background genes
dirname_read <- file.path(data_path, direction)
dn_aggr_signature_marray <- as.character(read_tsv(here("data/v2/signatures/microarray/aggregated_signatures/microarray_dn_aggregated_signature_0.8.tsv"))$GeneID)
enrichGO_res_dn_aggr_marray <- enrichGO(gene = dn_aggr_signature_marray,
                                        OrgDb = org.Hs.eg.db,
                                        readable = TRUE,
                                        ont = "BP",
                                        universe = bg_genes_combined_dn,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.1,
                                        minGSSize = 5,
                                        maxGSSize = 200)

enrichGO_res_dn_aggr_marray <- as.data.frame(enrichGO_res_dn_aggr_marray)

# Save results
output_file_dn_aggr_marray <- file.path(here("data/v2/microarray_pathways/dn"), paste0("GO_ORA_BGcorrected_dn_aggregated_microarray_0.8.tsv"))

write_tsv(enrichGO_res_dn_aggr_marray, file = output_file_dn_aggr_marray)
