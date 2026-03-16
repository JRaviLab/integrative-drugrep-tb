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

| Signature ID | Signature | Study | TB status | Origin type | Tissue origin | Cell/tissue type | Up genes | Down genes | Control samples | Disease samples |
|---|---|---|---|---|---|---|---|---|---|---|
| s1_d1 | GSE83456_GPL10558_EPTB_control_blood | GSE83456 | EPTB | primary | circulating | blood | 44 | 38 | 30 | 47 |
| s2_d1 | GSE83456_GPL10558_PTB_control_blood | GSE83456 | PTB | primary | circulating | blood | 54 | 91 | 31 | 45 |
| s3_d2 | GSE16250_GPL570_MTB_control_blood_pmbc | GSE16250 | MTB | primary | circulating | blood_PBMC | 56 | 13 | 3 | 3 |
| s4_d3 | GSE34151_GPL10558_MTB_control_blood_dcs | GSE34151 | MTB | primary | circulating | blood_dendritic_cells | 179 | 237 | 129 | 130 |
| s5_d4 | GSE139871_GPL10558_MTB_control_blood_monocyte | GSE139871 | MTB | primary | circulating | blood_monocyte | 11 | 32 | 4 | 20 |
| s6_d5 | GSE63548_GPL10558_LNTB_control_lymphnode | GSE63548 | LNTB | primary | other | lymphnode | 108 | 104 | 4 | 22 |
| s29_d24 | microarray_aggregated_signature | - | - | - | - | - | 98 | 98 | - | - |

---

### RNA-seq signatures

| Signature ID | Signature | Study | TB status | Origin type | Tissue origin | Cell/tissue type | Up genes | Down genes | Control samples | Disease samples |
|---|---|---|---|---|---|---|---|---|---|---|
| s7_d6 | GSE211974_MTB_control_blood_THP-1 | GSE211974 | MTB | cell line | circulating | THP-1 | 67 | 46 | 3 | 3 |
| s8_d7 | GSE141656_MTB_control_blood_blood | GSE141656 | MTB | primary | circulating | blood | 267 | 157 | 3 | 3 |
| s9_d8 | GSE193777_MTB_control_blood | GSE193777 | MTB | primary | circulating | blood | 215 | 108 | 36 | 10 |
| s10_d9 | GSE84076_MTB_control_blood | GSE84076 | MTB | primary | circulating | blood | 20 | 8 | 4 | 6 |
| s11_d10 | GSE148731_MTB_4h_blood_CD14plus_monocyte_derived_macrophages_MF1 | GSE148731 | MTB | primary | circulating | macrophage_MF1 | 50 | 7 | 3 | 3 |
| s12_d10 | GSE148731_MTB_24h_blood_CD14plus_monocyte_derived_macrophages_MF1 | GSE148731 | MTB | primary | circulating | macrophage_MF1 | 212 | 147 | 3 | 3 |
| s13_d10 | GSE148731_MTB_24h_blood_CD14plus_monocyte_derived_macrophages_MF2 | GSE148731 | MTB | primary | circulating | macrophage_MF2 | 286 | 219 | 4 | 3 |
| s14_d11 | GSE148171_MTB_control_blood_PBMC | GSE148171 | MTB | primary | circulating | PBMC | 176 | 52 | 4 | 5 |
| s15_d12 | GSE174566_MTB_control_blood_PBMC | GSE174566 | MTB | primary | circulating | PBMC | 209 | 95 | 18 | 18 |
| s16_d13 | GSE198557_MTB_control_blood_PBMC | GSE198557 | MTB | primary | circulating | PBMC | 156 | 48 | 6 | 6 |
| s17_d14 | GSE143627_MTB_control_blood_PBMC_macrophage | GSE143627 | MTB | primary | circulating | PBMC_macrophage | 102 | 83 | 3 | 3 |
| s18_d15 | GSE64179_MTB_control_blood_dendritic_cells | GSE64179 | MTB | primary | circulating | dendritic_cells | 271 | 243 | 6 | 6 |
| s19_d16 | GSE64182_MTB_control_blood_dendritic_cells | GSE64182 | MTB | primary | circulating | dendritic_cells | 235 | 88 | 3 | 3 |
| s20_d17 | GSE132283_MTB_control_2h_blood_macrophage | GSE132283 | MTB | primary | circulating | macrophage | 170 | 54 | 3 | 3 |
| s21_d17 | GSE132283_MTB_control_48h_blood_macrophage | GSE132283 | MTB | primary | circulating | macrophage | 336 | 243 | 3 | 3 |
| s22_d18 | GSE183912_MTB_control_blood_macrophage | GSE183912 | MTB | primary | circulating | macrophage | 25 | 2 | 4 | 4 |
| s23_d19 | GSE143731_MTB_control_blood_monocyte_derived_macrophages | GSE143731 | MTB | primary | circulating | monocyte_derived_macrophages | 6 | 7 | 4 | 4 |
| s24_d20 | GSE164287_MTB_control_blood_monocyte_derived_macrophages | GSE164287 | MTB | primary | circulating | monocyte_derived_macrophages | 81 | 84 | 4 | 4 |
| s25_d21 | GSE236156_MTB_control_blood_monocyte_derived_macrophage | GSE236156 | MTB | primary | circulating | monocyte_derived_macrophages | 213 | 161 | 11 | 6 |
| s26_d22 | GSE256184_MTB_control_0h_HEPG2 | GSE256184 | MTB | cell line | other | HepG2 | 167 | 124 | 4 | 4 |
| s27_d22 | GSE256184_MTB_48h_HEPG2 | GSE256184 | MTB | cell line | other | HepG2 | 410 | 277 | 4 | 4 |
| s28_d23 | GSE165708_MTB_control_BAL_bronchoalveolar_lavage | GSE165708 | MTB | primary | other | BAL_bronchoalveolar_lavage | 23 | 2 | 16 | 16 |
| s30_d24 | RNAseq_aggregated_signature | - | - | - | - | - | 98 | 98 | - | - |


## Drug data
We used drug perturbation data from the LINCS dataset GSE92742 level 5, implemented in the signatureSearchData resource and accessible through the signatureSearch R package.

## References

1. LINCS dataset GSE92742.  
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
2. signatureSearchData: Data package for signatureSearch.  
   https://bioconductor.org/packages/signatureSearchData/
3. signatureSearch: Environment for Gene Expression Searching and Functional Enrichment Analysis. 
   https://bioconductor.org/packages/signatureSearch/
