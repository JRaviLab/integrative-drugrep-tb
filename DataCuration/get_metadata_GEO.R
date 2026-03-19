# script to get metadata file given GEO accession ids
# created date: 07/13/23
# author: Kewalin Samart
# Modified by: Ling Thang
# last modified: 06/17/25
# This version processes studies with multiple platforms (GPLs)
# Usage: Rscript get_metadata_GEO.R "<GSE_IDs||PATHtoList>" <METADATA_DIR>

library(GEOquery)
library(readr)
library(dplyr) # Added for bind_rows()

## set up arguments
args <- commandArgs(TRUE)

# Check if args[1] is a file
if (file.exists(args[1])) {
  GSEs <- readLines(args[1]) # Reads from file
} else {
  GSEs <- strsplit(args[1], split = ",")[[1]] # Reads GSE IDs from a string
}

# args[2]: path to directory to save resulting metadata e.g., <SAVEPATH>
# If the path is relative, it will be relative to the current working directory.
dirname_meta <- args[2]

# Create the directory if it doesn't exist
if (!dir.exists(dirname_meta)) {
  dir.create(dirname_meta, recursive = TRUE)
}

# get dataset series matrix files
gse_list <- list()
for (id in GSEs) {
  print(paste("Downloading data for", id))
  gse_list[[id]] <- getGEO(id)
}

# obtain metadata for each platform and save to a single combined file per GSE

# Track success/failure for each GSE
success_gse <- c()
failed_gse <- c()

for (gse_id in names(gse_list)) {
  # Initialize a list to hold metadata from each platform for the current GSE
  all_platforms_metadata <- list()

  # gse_list[[gse_id]] is a list of ExpressionSet objects, one for each platform
  platforms <- gse_list[[gse_id]]
  for (eSet in platforms) {
    pd <- pData(eSet)

    # Add series_id (GSE) and platform_id (GPL) columns for traceability
    pd$series_id <- gse_id
    pd$platform_id <- annotation(eSet)

    # Add the platform's metadata dataframe to our list
    all_platforms_metadata <- append(
      all_platforms_metadata,
      list(as.data.frame(pd))
    )
  }

  # Combine all collected metadata dataframes into one.
  # bind_rows handles differing columns by filling missing values with NA.
  if (length(all_platforms_metadata) > 0) {
    combined_metadata <- dplyr::bind_rows(all_platforms_metadata)

    # Define the single output filename for the GSE
    output_filename <- file.path(
      dirname_meta,
      paste0("metadata_", gse_id, ".tsv")
    )

    print(paste("Saving combined metadata for", gse_id, "to:", output_filename))
    write_tsv(combined_metadata, output_filename)
    success_gse <- c(success_gse, gse_id)
  } else {
    print(paste("Warning: No metadata found for", gse_id))
    failed_gse <- c(failed_gse, gse_id)
  }
}

# Summary report
print("Script finished.")
if (length(failed_gse) > 0) {
  print(paste("Failed to retrieve:", paste(failed_gse, collapse = ", ")))
} else {
  print("All GSEs were successfully retrieved.")
}
