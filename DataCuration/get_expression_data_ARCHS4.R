# script to download selected samples
# Copy code and run on a local machine to initiate download
# Original Author : Kawalin Samart
# Modified by: Ling Thang
# Date : Oct 14 2025
# Usage: Rscript get_expression_data_ARCHS4.R <GSE_ids_file> <output_dir>

library("rhdf5") # can be installed using Bioconductor
library(readr)

# Function to process each GSE ID
process_gse <- function(gse_id, destination_file, output_dir) {
  # expression file path
  extracted_expression_file <- paste0(
    output_dir,
    "/", gse_id, ".tsv"
  )

  # Get disease name from output_dir (assumes format "data/expressions/DISEASE")
  disease_name <- basename(output_dir)

  # Construct metadata file path dynamically
  sampids_path <- file.path(
    "data/metadata",
    disease_name,
    paste0("metadata_", gse_id, ".tsv")
  )

  sampids <- read.csv(sampids_path, sep = "\t", row.names = NULL)

  # Selected samples to be extracted
  samp <- sampids$geo_accession

  # Retrieve information from compressed data
  samples <- h5read(destination_file, "meta/samples/geo_accession")
  genes <- h5read(destination_file, "meta/genes/symbol")

  # Identify columns to be extracted
  sample_locations <- which(samples %in% samp)

  # Extract gene expression from compressed data
  expression <- t(h5read(
    destination_file, "data/expression",
    index = list(sample_locations, seq_along(genes))
  ))
  H5close()
  rownames(expression) <- genes
  colnames(expression) <- samples[sample_locations]

  # Try writing the file with error handling
  tryCatch(
    {
      write.table(expression,
        file = extracted_expression_file,
        sep = "\t", quote = FALSE, col.names = NA
      )

      # Check if file exists before confirming success
      if (file.exists(extracted_expression_file)) {
        cat(paste(gse_id, "Successfully processed!\n"))
      } else {
        cat(paste("WARNING: File was not created!", gse_id))
      }
    },
    error = function(e) {
      cat(paste("ERROR: Failed to write file:", gse_id))
      cat(paste("Error message:", e$message))
    }
  )
}

# Main script
args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("Please provide the path to the GSE IDs file and the output directory.")
}

gse_ids_file <- args[1]
output_dir <- args[2]
destination_file <- "data/human_gene_v2.5.h5"

# Read GSE IDs from the specified file from tsv
gse_ids <- unique(read_tsv(gse_ids_file, col_names = FALSE)[[1]])

# Process each GSE -- with status tracking
all_success <- TRUE

for (gse_id in gse_ids) {
  tryCatch(
    {
      process_gse(gse_id, destination_file, output_dir)
    },
    error = function(e) {
      cat(paste("ERROR processing", gse_id, ":", e$message))
      all_success <<- FALSE # Mark failure
    }
  )
}

# Print final completion message based on success
if (all_success) {
  cat("All GSE IDs have been successfully processed!")
} else {
  cat("Completed with errors. Check logs for details.")
}
