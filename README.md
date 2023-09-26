# drugrep_baseline_comparison
Baseline comparison for drug repurposing project
## Databases
- [LINCS](clue.io) –– drug expression data; access/obtain using [`signatureSearch`](https://github.com/girke-lab/signatureSearch) R package 
- [GEO](https://www.ncbi.nlm.nih.gov/geo/) –– disease expression data; access/obtain using [`ARCHS4py`](https://github.com/MaayanLab/archs4py) for [RNAseq data](https://maayanlab.cloud/archs4/) and [`refine.bio`](https://www.refine.bio/) for microarray data

## Steps
1. Getting drug control samples
   - use `signatureSearch` to look into drug cell line information then record the names of cell lines with cell_type == 'normal'
   - extract the expression data for the control samples obtained from the previous step i.e., all the available untreated cell line in LINCS and CMAP drug databases
3. Cleaning up data
   - make sure to use `Gene ID` i.e., `Ensembl` throughout; using mapping file in the `annotaion` folder
   - for both disease and drug data, subset the expression data to include only "landmark genes" –– 978 genes total
   - make sure that the expression data for all the samples have the same order of genes (uniform row names) 
4. For all pairs of disease- and drug-control samples, calculate correlation starting with pearson, spearman etc. So, in the end, we would have a correlation matrix for each similarity matrix where rows and columns are disease- and drug-control samples (or the other way around). 
