# Baseline analysis
This section describes how to reproduce the baseline analysis used in this study.

## Overview

The baseline analysis compares the disease control samples with drug control samples from the LINCS L1000 dataset. We use LINCS Level 3 expression data, which contains normalized and inferred gene expression values.
The analysis consists of three steps:

1. Download LINCS Level 3 expression data and metadata
2. Extract drug control samples and restrict to the L1000 gene space
3. Compare disease and drug control samples 

## Requirements (Python packages)

The pipeline requires:

- python >= 3.8
- pandas
- numpy
- cmapPy
- h5py

Install dependencies with:
```
pip install cmapPy pandas numpy h5py
```

### Data sources

Drug perturbation data are obtained from the LINCS L1000 dataset.

Dataset: GEO accession GSE92742

Documentation:

- https://clue.io/connectopedia/guide_to_geo_l1000_data
- https://clue.io/connectopedia/gctx_format

## 1. Retrieve LINCS level 3 expression drug data and metadata
### LINCS level 3 expression data
  - Access/obtain using [`cmapPy`](https://github.com/cmap/cmapPy) Python package; [descriptive details](https://clue.io/connectopedia/guide_to_geo_l1000_data) on LINCS data.
  - Since we are extracting control samples only, we will be using LINCS data level 3 which is normalized and inferred. The data could be downloaded from [this GEO site]([https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)) or using the command lines below.

  ```
  mkdir -p data/baseline_analysis && \
  wget -P data/baseline_analysis \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz
  ```
#### Notes

- This file is ~50 GB compressed, so downloading may take significant time and disk space.
- Using an HPC environment is recommended for faster processing.
- The `wget -c` flag allows interrupted downloads to resume.

### LINCS level 3 metadata
  ```
  mkdir -p data/metadata && \
  wget -O data/metadata/GSE92742_Broad_LINCS_inst_info.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_inst_info.txt.gz
  ```
  This metadata includes:
  - perturbation type
  - compound identifiers
  - cell line
  - dose
  - time point
  - instance IDs

These metadata are required to identify control samples and drug perturbation experiments.

## 2. Extract drug control samples
Run the following script to extract control samples from LINCS experiments and restrict the data to the L1000 gene space.
```
python scripts/suppl_baseline_get_drugdata_lincslevel3.py
```

This script performs:

- loading LINCS Level 3 expression data
- filtering for control samples
- mapping metadata to expression profiles
- restricting expression to L1000 genes
- saving processed drug control data for downstream analysis

## 3. Baseline comparison analysis
To reproduce the baseline comparison and Figure S9, run the notebook:
```
figures/figureS9/suppl_baseline_comparison.ipynb
```

This notebook performs:

- comparison between disease control samples and drug control samples
- baseline similarity analysis
- visualization of results used to generate Figure S9
