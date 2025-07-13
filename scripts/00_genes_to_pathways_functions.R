# script to perform pathway analyses: GO Over-represenation analysis
# last modified: 07/13/25
# Kewalin Samart

library(tidyverse)
library(here)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(org.Hs.eg.db)

source(here("scripts/00_background_genes_PA_functions.R"))

##### PARAMETERS #####
technologies <- c("microarray", "RNAseq")
directions <- c("up", "dn")

##### FUNCTION: Split DE Results #####
split_DE_results <- function(technology) {

  metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
  signature_info <- read_tsv(metadata_path) %>% filter(signature == 1)

  data_path <- here(paste0("data/v2/DE_results/", technology))
  dir.create(file.path(data_path, "up"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(data_path, "dn"), showWarnings = FALSE, recursive = TRUE)

  filenames <- list.files(data_path, pattern = "\\.tsv$", full.names = FALSE)
  signature_filenames <- filenames[grepl(paste(signature_info$series_id, collapse = "|"), filenames)]

  for (file in signature_filenames) {
    de_df <- read_tsv(file.path(data_path, file))

    up_genes <- de_df %>% filter(adj.P.Val < 0.05, log2FoldChange > 0)
    dn_genes <- de_df %>% filter(adj.P.Val < 0.05, log2FoldChange < 0)

    base_filename <- sub("\\.tsv$", "", file)
    write_tsv(up_genes, file.path(data_path, "up", paste0(base_filename, ".up.tsv")))
    write_tsv(dn_genes, file.path(data_path, "dn", paste0(base_filename, ".dn.tsv")))
  }
}


##### FUNCTION: Pathway Enrichment #####
run_pathway_analysis <- function(technology, direction) {

  metadata_path <- here(paste0("data/v2/signatures/", technology, "_TB_signature_run_info.tsv"))
  data_path <- here(paste0("data/v2/DE_results/", technology))
  output_dir <- here(paste0("data/v2/", technology, "_pathways/", direction))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  input_dir <- file.path(data_path, direction)
  filenames <- list.files(input_dir, pattern = paste0("\\.", direction, "\\.tsv$"))

  # Background genes
  extension <- paste0("_", direction, ".tsv")
  bg_genes_raw <- get_bg_genes(bg_source = "input data",
                               metadata_path = metadata_path,
                               data_path = data_path,
                               extension = extension)

  bg_genes <- intersect(bg_genes_raw, GO_genes())

  # Iterate over DE result files
  for (file in filenames) {
    genes_df <- read_tsv(file.path(input_dir, file))

    gene_list <- as.character(genes_df$GeneID)

    enrichGO_res <- enrichGO(gene = gene_list,
                             OrgDb = org.Hs.eg.db,
                             readable = TRUE,
                             ont = "BP",
                             universe = bg_genes,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.1,
                             minGSSize = 5,
                             maxGSSize = 200) %>%
      as.data.frame()

    base_filename <- sub(paste0("\\.", direction, "\\.tsv$"), "", file)
    output_file <- file.path(output_dir, paste0("GO_ORA_BGcorrected_", base_filename, ".tsv"))

    write_tsv(enrichGO_res, file = output_file)
  }
}


##### PIPELINE #####
for (technology in technologies) {
  message("Splitting DE results for: ", technology)
  split_DE_results(technology)

  for (direction in directions) {
    message("Running pathway analysis for: ", technology, " (", direction, ")")
    run_pathway_analysis(technology, direction)
  }
}

