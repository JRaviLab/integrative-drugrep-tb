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
Details of disease data can be found in `DataCuration/`.  
The summary of disease signatures used in this study is provided in [**Table S1**](https://github.com/JRaviLab/integrative-drugrep-tb/blob/main/figures/tableS1/tableS1_signature_summary_table.pdf).

### Microarray signatures

| No. | ID | Signature name | TB status | Tissue origin | Cell type | Origin type | Dataset type | Upregulated genes | Downregulated genes | Control samples | Disease samples |
|-----|----|---------------|----------|--------------|----------|-------------|--------------|------------------|--------------------|----------------|----------------|
| 1 | s1_d1 | GSE83456_GPL10558_s1_d1 | PTB | Blood | Whole blood | Primary | Circulating | 54 | 91 | 31 | 45 |
| 2 | s2_d1 | GSE83456_GPL10558_s2_d1 | EPTB | Blood | Whole blood | Primary | Circulating | 44 | 38 | 30 | 47 |
| 3 | s3_d2 | GSE16250_GPL570_s3_d2 | MTB | Blood | Mononuclear cell | Primary | Circulating | 56 | 13 | 3 | 3 |
| 4 | s4_d3 | GSE34151_GPL10558_s4_d3 | MTB | Blood | Dendritic cell | Primary | Circulating | 179 | 237 | 129 | 130 |
| 5 | s5_d4 | GSE139871_GPL10558_s5_d4 | MTB | Blood | Monocyte | Primary | Circulating | 11 | 32 | 4 | 20 |
| 6 | s6_d5 | GSE63548_GPL10558_s6_d5 | LNTB | Lymph node | Lymph node | Primary | Other | 108 | 104 | 4 | 22 |
| 29 | s29_d24 | microarray_aggregated_signature | - | - | - | - | - | 98 | 98 | - | - |

---

### RNA-seq signatures

| No. | ID | Signature name | TB status | Tissue origin | Cell type | Origin type | Dataset type | Upregulated genes | Downregulated genes | Control samples | Disease samples |
|-----|----|---------------|----------|--------------|----------|-------------|--------------|------------------|--------------------|----------------|----------------|
| 7 | s7_d6 | GSE211974_s7_d6 | MTB | THP1 | Macrophage | Cell line | Circulating | 67 | 46 | 3 | 3 |
| 8 | s8_d7 | GSE141656_s8_d7 | MTB | Blood | Whole blood | Primary | Circulating | 267 | 157 | 3 | 3 |
| 9 | s9_d8 | GSE193777_s9_d8 | MTB | Blood | Whole blood | Primary | Circulating | 215 | 108 | 36 | 10 |
| 10 | s10_d9 | GSE84076_s10_d9 | MTB | Blood | Whole blood | Primary | Circulating | 20 | 8 | 4 | 6 |
| 11 | s11_d10 | GSE148731_s11_d10 | MTB | Blood | M1 macrophage | Primary | Circulating | 212 | 147 | 3 | 3 |
| 12 | s12_d10 | GSE148731_s12_d10 | MTB | Blood | M1 macrophage | Primary | Circulating | 50 | 7 | 3 | 3 |
| 13 | s13_d10 | GSE148731_s13_d10 | MTB | Blood | M2 macrophage | Primary | Circulating | 286 | 219 | 4 | 3 |
| 14 | s14_d11 | GSE148171_s14_d11 | MTB | Blood | Mononuclear cell | Primary | Circulating | 176 | 52 | 4 | 5 |
| 15 | s15_d12 | GSE174566_s15_d12 | MTxjB | Blood | PBMC | Primary | Circulating | 209 | 95 | 18 | 18 |
| 16 | s16_d13 | GSE198557_s16_d13 | MTB | Blood | PBMC | Primary | Circulating | 156 | 48 | 6 | 6 |
| 17 | s17_d14 | GSE143627_s17_d14 | MTB | Blood | Mononuclear leukocyte | Primary | Circulating | 102 | 83 | 3 | 3 |
| 18 | s18_d15 | GSE64179_s18_d15 | MTB | Blood | Dendritic cell | Primary | Circulating | 271 | 243 | 6 | 6 |
| 19 | s19_d16 | GSE64182_s19_d16 | MTB | Blood | Dendritic cell | Primary | Circulating | 235 | 88 | 3 | 3 |
| 20 | s20_d17 | GSE132283_s20_d17 | MTB | Blood | Macrophage | Primary | Circulating | 170 | 54 | 3 | 3 |
| 21 | s21_d17 | GSE132283_s21_d17 | MTB | Blood | Macrophage | Primary | Circulating | 336 | 243 | 3 | 3 |
| 22 | s22_d18 | GSE183912_s22_d18 | MTB | Blood | Macrophage | Primary | Circulating | 25 | 2 | 4 | 4 |
| 23 | s23_d19 | GSE143731_s23_d19 | MTB | Blood | Macrophage | Primary | Circulating | 6 | 7 | 4 | 4 |
| 24 | s24_d20 | GSE164287_s24_d20 | MTB | Blood | Macrophage | Primary | Circulating | 81 | 84 | 4 | 4 |
| 25 | s25_d21 | GSE236156_s25_d21 | MTB | Blood | Macrophage | Primary | Circulating | 213 | 161 | 11 | 6 |
| 26 | s26_d22 | GSE256184_s26_d22 | MTB | HepG2 | Hepatocyte | Cell line | Other | 410 | 277 | 4 | 4 |
| 27 | s27_d22 | GSE256184_s27_d22 | MTB | HepG2 | Hepatocyte | Cell line | Other | 167 | 124 | 4 | 4 |
| 28 | s28_d23 | GSE165708_s28_d23 | MTB | BALF | Alveolar system | Primary | Other | 23 | 2 | 16 | 16 |
| 30 | s30_d24 | RNAseq_aggregated_signature | - | - | - | - | - | 98 | 98 | - | - |


## Drug data
We used drug perturbation data from the LINCS dataset GSE92742 level 5, implemented in the signatureSearchData resource and accessible through the signatureSearch R package.

## References

1. LINCS dataset GSE92742.  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
2. signatureSearchData: Data package for signatureSearch.  
   https://bioconductor.org/packages/signatureSearchData/
3. signatureSearch: Environment for Gene Expression Searching and Functional Enrichment Analysis. 
   https://bioconductor.org/packages/signatureSearch/
