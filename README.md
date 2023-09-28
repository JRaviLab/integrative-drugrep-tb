# drugrep_baseline_comparison
Baseline comparison for drug repurposing project
## Databases
### Drug data
- [LINCS](clue.io) –– drug expression data; access/obtain using [`signatureSearch`](https://github.com/girke-lab/signatureSearch) R package; [descriptive details](https://clue.io/connectopedia/guide_to_geo_l1000_data) on LINCS data. Since we are extracting control samples only, we will be using LINCS data level 3 which is normalized and inferred. The data could be manually downloaded from [this GEO site](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) or using the command line below
  
  **[GCTX format](https://clue.io/connectopedia/gctx_format) (expression data + metadata)**
  ```
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz
  ```
  ![GCTx example](https://clue.io/connectopedia/images/gctx_format_images/image_0.png)

  **[GCT format](https://clue.io/connectopedia/gct_format) (expression data; parse(GCTX) ––> GCT object; GCT@mat ––> expression matrix )**
  ```
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n113012x22268_2015-12-31.gct.gz
  ```
  ![GCT example](https://clue.io/connectopedia/images/gct_format_images/image_0.png)

  **LINCS data level 3 metadata**

  Here we get `inst_id` corresponding to `pert_type` == "ctl_vehicle" (control samples)
  ```
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz
  ```

  Example data:
  ```
  inst_id	                      cell_id	  det_plate	               det_well	  pert_mfc_id	  pert_dose	  pert_dose_unit	pert_id	  pert_iname	pert_type	  pert_time	pert_time_unit
  LJP005_A375_24H_X1_B19:A03	  A375	    LJP005_A375_24H_X1_B19	 A03	      DMSO	        -666.0	    -666	          DMSO	    DMSO	      ctl_vehicle	24.0      h
  LPROT001_A375_6H_X1_B20:B03	  A375	    LPROT001_A375_6H_X1_B20	 B03	      DMSO	        -666.0	    -666	          DMSO	    DMSO	      ctl_vehicle	6.0       h
  LPROT001_A375_6H_X1_B20:B05	  A375	    LPROT001_A375_6H_X1_B20	 B05	      DMSO	        -666.0	    -666	          DMSO	    DMSO	      ctl_vehicle	6.0       h
  LPROT002_A375_6H_X1_B22:B03	  A375	    LPROT002_A375_6H_X1_B22	 B03	      DMSO	        -666.0	    -666	          DMSO	    DMSO	      ctl_vehicle	6.0       h
  ```

### Disease data
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
