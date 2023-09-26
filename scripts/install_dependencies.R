# run this script to install R dependencies for the project
# to be modified as needed
# created date: 09/26/23
# Kewalin Samart

# required R packages from CRAN: https://cran.r-project.org/
install.packages("stringr", "tidyverse", "purrr")

# required R packages from Bioconductor: https://bioconductor.org/
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("biomaRt",
                       "limma",
                       "oligo",
                       "GEOquery",
                       "GEOmetadb",
                       "affycoretools",
                       "biomaRt",
                       "qvalue",
                       "signatureSearch", # extract untreated drug expression data
                       "Rcpp", 
                       "SummarizedExperiment", 
                       "org.Hs.eg.db"))