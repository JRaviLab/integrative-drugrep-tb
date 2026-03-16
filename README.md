# Integrative Transcriptome-Based Drug Repurposing in Tuberculosis

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
- Python 3 (for `scripts/00_cleanup_expression_data.ipynb`)

### Setup

```r
# Install renv if needed
install.packages("renv")

# Restore the R package library
renv::restore()
```

### Running the pipeline

Scripts are designed to be run sequentially (steps 00 → 05). Each main script accepts command-line arguments; defaults are provided for the RNAseq/LINCS analysis path. Example:

```bash
Rscript scripts/02_drugrep_get_prediction.R \
  data/signatures/RNASeq_TB_signature_run_info.tsv \
  data/signatures/RNAseq \
  LINCS LINCS \
  results/RNAseq/LINCS
```

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
