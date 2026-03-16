# Data
This directory contains processed datasets and intermediate results used in the analysis pipeline for the study.

The folder structure is organized to separate expression datasets, differential expression results, pathway analyses, and metadata required for reproducing the analyses. 
```
.
├── DE_results/                # Differential expression results
├── RNAseq_data_forDE/         # Processed RNA-seq datasets used for DE analysis
├── RNAseq_pathways/           # Pathway analysis results for RNA-seq datasets
├── baseline_analysis/         # Baseline drug-control comparison analysis
├── metadata/                  # Dataset metadata and annotations
├── microarray_data_forDE/     # Processed microarray datasets used for DE analysis
├── microarray_pathways/       # Pathway analysis results for microarray datasets
├── signatures/                # Final disease signatures used in analysis
└── README.md
```
## Disease data
Details of disease data can be found in `DataCuration/.`

## Drug data
We used drug perturbation data from the LINCS dataset GSE92742 level 5, implemented in the signatureSearchData resource and accessible through the signatureSearch R package.

## References

1. LINCS dataset GSE92742.  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
2. signatureSearchData: Data package for signatureSearch.  
   https://bioconductor.org/packages/signatureSearchData/
3. signatureSearch: Environment for Gene Expression Searching and Functional Enrichment Analysis. 
   https://bioconductor.org/packages/signatureSearch/
