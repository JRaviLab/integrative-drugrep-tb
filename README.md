# Integrative transcriptome-based drug repurposing in tuberculosis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the analysis code and processed data for the study:

> Samart K, Thang L, Buskirk L, Tonielli A, Krishnan A\*, Ravi J\*. *Integrative transcriptome-based drug repurposing in tuberculosis.* bioRxiv (2025). https://doi.org/10.1101/2025.06.02.657296
>
> \* co-corresponding authors

We integrate transcriptomic signatures from multiple TB microarray and RNA-seq datasets with drug–gene connectivity databases (CMap, LINCS) to identify and prioritize drug repurposing candidates for tuberculosis.

---

## Repository structure

```
integrative-drugrep-tb/
├── scripts/        # R and Python analysis scripts (numbered by step)
├── vignette/       # Quarto documents with worked examples for each step
├── figures/        # Code and outputs for manuscript figures (figure1–5, figureS1–S9)
├── data/           # Input data: DE results, signatures, and metadata
├── results/        # Pipeline outputs: connectivity scores and drug rankings
└── renv.lock       # R package lockfile for reproducibility
```

---

## Workflow

The analysis follows a five-step pipeline:

| Step | Scripts | Description |
|------|---------|-------------|
| **00** | `00_multids_microarray_DEwithlimma.R`, `00_multids_RNAseq_DEwithDESeq2.R`, `00_cleanup_expression_data.ipynb` | Differential expression analysis across TB datasets |
| **01** | `01_signature_aggregation_functions.R`, `01_signature_landmark_prep_functions.R` | Aggregate DE signatures; filter to L1000 landmark genes |
| **02** | `02_drugrep_get_prediction.R`, `02_signatureSearch_connectivity_scores_functions.R` | Score drug candidates via CMap 1.0, LINCS, and correlation methods |
| **03** | `03_summarize_drugs_methodswise_functions.R` | Summarize and compare predictions across scoring methods |
| **04** | `04_*_RankAggregation_*.R` | Aggregate drug rankings across methods and datasets |
| **05** | `05_high_confidence_drug_prediction.R`, `05_build_drugtarget_network_functions.R` | Identify high-confidence candidates and build drug–target networks |

See [`scripts/README.md`](scripts/README.md) for a full table of scripts and their purposes, and [`vignette/`](vignette/) for step-by-step worked examples.

---

## Reproducing the analysis

### Requirements

- R 4.4.2 with Bioconductor 3.20
- Python 3 (for [baseline analysis](https://github.com/JRaviLab/integrative-drugrep-tb/tree/main/figures/figureS9))

### Setup

```r
# Install renv if needed
install.packages("renv")

# Restore the R package library
renv::restore()
```

### Running the pipeline
### 1. Generate disease signatures via differential expression analysis
#### 1.1 Microarray
```
Rscript scripts/00_multids_microarray_DEwithlimma.R <metadata_file.tsv> <padj_cutoff>
```
Arguments:
- `meta_class_file.tsv`  : Tab‑separated file listing all datasets: `./data/microarray_data_forDE/clean_TB_sample_metadata_classification.tsv`.

  Mandatory columns:
  - `series_id`       (GEO study identifier)
  - `geo_accession`   (GEO sample identifier )
  - `SIGNATURE_NAME`  (name of signature containing unique sample conditions)
  - `EXPRMAT_PATH`    (path to raw‑count matrix TSV)
  - `CLASSIFICATION`  (labels:  `disease_without_treatment`  or  `healthy_control_without_treatment`)
- `padj_cutoff`    : Adjusted‑p significance threshold (default 0.05)

#### 1.2 RNAseq
```
Rscript scripts/00_multids_RNAseq_DESeq2.R <metadata_file.tsv> <padj_cutoff>
```
Arguments:
- `meta_class_file.tsv`  : Tab‑separated file listing all datasets: `./data/RNAseq_data_forDE/clean_TB_sample_metadata_classification.tsv`.

  Mandatory columns:
  - `series_id`       (GEO study identifier)
  - `geo_accession`   (GEO sample identifier )
  - `SIGNATURE_NAME`  (name of signature containing unique sample conditions)
  - `EXPRMAT_PATH`    (path to raw‑count matrix TSV)
  - `CLASSIFICATION`  (labels:  `disease_without_treatment`  or  `healthy_control_without_treatment`)
- `padj_cutoff`    : Adjusted‑p significance threshold (default 0.05)

### 2. Prioritize drug candidates using multiple connectivity scores
Example run for quantifying candidate drugs reversing RNAseq TB signatures using CMAP 2.0 methods (i.e., LINCS):
```bash
Rscript scripts/02_drugrep_get_prediction.R \
  data/signatures/RNASeq_TB_signature_run_info.tsv \
  data/signatures/RNAseq \
  LINCS LINCS \
  results/RNAseq/LINCS
```
Arguments:
- `sig_metadata_path` – Path to the RNA-seq signature metadata file (default: data/signatures/RNASeq_TB_signature_run_info.tsv).
- `sig_data_path` – Directory containing RNA-seq signature data files (default: data/signatures/RNAseq).
- `drugdb_name` – Drug perturbation database to use. Options: LINCS, CMAP (default: LINCS).
- `score_method` – Method used to compute signature similarity scores. Options: LINCS, CMAP, Cor_spearman, Cor_pearson (default: LINCS).
- `output_dir` – Directory where output results will be saved (default: results/RNAseq/LINCS).


### 3. Summarize drug predictions



### 4.
Additional downstream analyses and figure generation can be reproduced using the Quarto (.qmd) notebooks provided in the `vignettes/` directory.

---

## Citation

If you use this code or data, please cite:

```
Samart K, Thang L, Buskirk L, Tonielli A, Krishnan A*, Ravi J*.
Integrative transcriptome-based drug repurposing in tuberculosis.
bioRxiv (2025). https://doi.org/10.1101/2025.06.02.657296
(* co-corresponding authors)
```

A `CITATION.cff` file is also provided for automated citation tools.

---

## License

MIT © 2026 [JRaviLab](https://jravilab.github.io/) & [KrishnanLab](https://www.thekrishnanlab.org). See [LICENSE](LICENSE).
