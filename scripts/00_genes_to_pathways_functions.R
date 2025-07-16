# script to perform pathway analyses: GO/KEGG/Reactome ORA
# last modified: 07/13/25
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

technologies <- c("microarray", "RNAseq")
directions <- c("up", "dn")

for (technology in technologies) {

  metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
  data_path <- here(paste0("data/v2/DE_results/",technology))

  for (direction in directions) {

    # Set output directory
    dirname <- here(paste0("data/v2/", technology, "_pathways/", direction))
    if (!dir.exists(dirname)) {
      dir.create(dirname, recursive = TRUE)
    }

    # Read gene files from up/dn folder
    dirname_read <- file.path(data_path, direction)
    filenames <- list.files(path = dirname_read, pattern = paste0(direction, ".tsv"))

    # Define background genes
    if(technology == "RNAseq") {
      extension <- paste0("_DESeq2_", direction, ".tsv")
    } else if(technology == "microarray") {
      extension <- paste0("_limma_", direction, ".tsv")
    }

    extra_arg <- NULL  # This seems unused

    bg_genes_source <- "input data"
    bg_genes_ <- get_bg_genes(bg_source = bg_genes_source,
                              metadata_path = metadata_path,
                              data_path = dirname_read,
                              extension = extension)

    bg_genes <- intersect(bg_genes_, GO_genes())  # Restrict to GO-annotated genes

    # Iterate through each file
    itr <- 1
    for (file in filenames) {

      print(paste0("Running ", technology, " ", direction, " iteration: ", itr))
      itr <- itr + 1

      genes_df <- read_tsv(file.path(dirname_read, file))

      # Use gene IDs as input for ORA
      gene_vector <- genes_df$log2FoldChange
      names(gene_vector) <- as.character(genes_df$GeneID)

      gene_list <- names(gene_vector)

      # Run ORA
      print(paste0("Getting enriched GO pathways for ", file))

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
bg_genes_up_rnaseq <- get_bg_genes(bg_source = bg_genes_source,
                          metadata_path = metadata_path,
                          data_path = dirname_read,
                          extension = extension)
bg_genes_up_rnaseq <- intersect(bg_genes_up_rnaseq, GO_genes())  # Restrict to GO-annotated genes
up_aggr_signature_rnaseq <- as.character(read_tsv(here("data/v2/signatures/RNAseq/aggregated_signatures/up_aggregated_signature.tsv"))$GeneID)
enrichGO_res_up_aggr_rnaseq <- enrichGO(gene = up_aggr_signature_rnaseq,
                         OrgDb = org.Hs.eg.db,
                         readable = TRUE,
                         ont = "BP",
                         universe = bg_genes_up_rnaseq,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         minGSSize = 5,
                         maxGSSize = 200)

enrichGO_res_up_aggr_rnaseq <- as.data.frame(enrichGO_res_up_aggr_rnaseq)

# Save results
output_file_up_aggr_rnaseq <- file.path(here("data/v2/RNAseq_pathways/up"), paste0("GO_ORA_BGcorrected_up_aggregated_RNAseq.tsv"))

write_tsv(enrichGO_res_up_aggr_rnaseq, file = output_file_up_aggr_rnaseq)


### dn signatures
direction <- "dn"
if(technology == "RNAseq"){extension <- paste0("_DESeq2_",direction,".tsv")}
if(technology == "microarray"){extension <- paste0("_limma_",direction,".tsv")}
# Define background genes
dirname_read <- file.path(data_path, direction)
bg_genes_dn_rnaseq <- get_bg_genes(bg_source = bg_genes_source,
                            metadata_path = metadata_path,
                            data_path = dirname_read,
                            extension = extension)
bg_genes_dn_rnaseq <- intersect(bg_genes_dn_rnaseq, GO_genes())  # Restrict to GO-annotated genes
dn_aggr_signature_rnaseq <- as.character(read_tsv(here("data/v2/signatures/RNAseq/aggregated_signatures/dn_aggregated_signature.tsv"))$GeneID)
enrichGO_res_dn_aggr_rnaseq <- enrichGO(gene = dn_aggr_signature_rnaseq,
                                        OrgDb = org.Hs.eg.db,
                                        readable = TRUE,
                                        ont = "BP",
                                        universe = bg_genes_dn_rnaseq,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.1,
                                        minGSSize = 5,
                                        maxGSSize = 200)

enrichGO_res_dn_aggr_rnaseq <- as.data.frame(enrichGO_res_dn_aggr_rnaseq)

# Save results
output_file_dn_aggr_rnaseq <- file.path(here("data/v2/RNAseq_pathways/dn"), paste0("GO_ORA_BGcorrected_dn_aggregated_RNAseq.tsv"))

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
bg_genes_up_marray <- get_bg_genes(bg_source = bg_genes_source,
                                   metadata_path = metadata_path,
                                   data_path = dirname_read,
                                   extension = extension)
bg_genes_up_marray <- intersect(bg_genes_up_marray, GO_genes())  # Restrict to GO-annotated genes
up_aggr_signature_marray <- as.character(read_tsv(here("data/v2/signatures/microarray/aggregated_signatures/up_aggregated_signature.tsv"))$GeneID)
enrichGO_res_up_aggr_marray <- enrichGO(gene = up_aggr_signature_marray,
                                        OrgDb = org.Hs.eg.db,
                                        readable = TRUE,
                                        ont = "BP",
                                        universe = bg_genes_up_marray,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.1,
                                        minGSSize = 5,
                                        maxGSSize = 200)

enrichGO_res_up_aggr_marray <- as.data.frame(enrichGO_res_up_aggr_marray)

# Save results
output_file_up_aggr_marray <- file.path(here("data/v2/microarray_pathways/up"), paste0("GO_ORA_BGcorrected_up_aggregated_microarray.tsv"))

write_tsv(enrichGO_res_up_aggr_marray, file = output_file_up_aggr_marray)

### dn signatures
direction <- "dn"
if(technology == "RNAseq"){extension <- paste0("_DESeq2_",direction,".tsv")}
if(technology == "microarray"){extension <- paste0("_limma_",direction,".tsv")}
# Define background genes
dirname_read <- file.path(data_path, direction)
bg_genes_dn_marray <- get_bg_genes(bg_source = bg_genes_source,
                                   metadata_path = metadata_path,
                                   data_path = dirname_read,
                                   extension = extension)
bg_genes_dn_marray <- intersect(bg_genes_dn_marray, GO_genes())  # Restrict to GO-annotated genes
dn_aggr_signature_marray <- as.character(read_tsv(here("data/v2/signatures/microarray/aggregated_signatures/dn_aggregated_signature.tsv"))$GeneID)
enrichGO_res_dn_aggr_marray <- enrichGO(gene = dn_aggr_signature_marray,
                                        OrgDb = org.Hs.eg.db,
                                        readable = TRUE,
                                        ont = "BP",
                                        universe = bg_genes_dn_marray,
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.1,
                                        minGSSize = 5,
                                        maxGSSize = 200)

enrichGO_res_dn_aggr_marray <- as.data.frame(enrichGO_res_dn_aggr_marray)

# Save results
output_file_dn_aggr_marray <- file.path(here("data/v2/microarray_pathways/dn"), paste0("GO_ORA_BGcorrected_dn_aggregated_microarray.tsv"))

write_tsv(enrichGO_res_dn_aggr_marray, file = output_file_dn_aggr_marray)
