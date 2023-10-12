# drugrep_baseline_comparison
Baseline comparisons for drug repurposing project
## Databases
### Drug data
- [LINCS](clue.io) –– drug expression data
  - Access/obtain using [`cmapPy`](https://github.com/cmap/cmapPy) Python package; [descriptive details](https://clue.io/connectopedia/guide_to_geo_l1000_data) on LINCS data.
  - Since we are extracting control samples only, we will be using LINCS data level 3 which is normalized and inferred. The data could be manually downloaded from [this GEO site](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) or using the command lines below.
  
  **[GCTX format](https://clue.io/connectopedia/gctx_format) (expression data + metadata)**
  ```
  wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz
  ```
  ![GCTx example](https://clue.io/connectopedia/images/gctx_format_images/image_0.png)

  **[GCT format](https://clue.io/connectopedia/gct_format) (expression data; parse(GCTx) ––> GCT object; GCT@mat ––> expression matrix )**
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
   Drug data:
   - Data wrangling: look into the metadata file: `GSE70138_Broad_LINCS_inst_info_2017-03-06.txt` then record the `inst_id` with `pert_type` == 'normal'
   - Get all the available untreated cell line in LINCS and CMAP drug databases: get drug control samples: use `cmapPy` to extract the expression data of the control samples obtained from the previous step from the gct object (parsed GCTx file: `GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz`).
   Disease data:
  - Only use control samples. To identify which samples are control samples, use the matadata where we have sample description to refer.
3. Cleaning up data
   - Make sure to use `Gene ID` i.e., `Entrezid` throughout; using mapping file in the `annotaion` folder
   - For both disease and drug data, subset the expression data to include only "landmark genes" –– 978 genes total
   - Dealing with duplicates: (i) take the one with maximum expression value, (ii) use the mean expression of the gene, or (iii) use the median expression of the gene. -- do this per sample.
4. Data exploration
   - Record size of intersecting genes for each pairwise of disease-drug control samples (could be a heatmap...). ** pay attention to the band with no overlaps; go back and check which datasets (GSE..) the samples are coming from.
   - Plot a histrogram of the overlapping sizes of gene sets (x-axis: intersecting size, y-axis: counts; should look very left skewed meaning the genes in most of the sample pairs highly overlap). 
5. Result format
   - Correlation matrices: For all pairs of disease- and drug-control samples, calculate correlation matrics starting with Pearson, Spearman, and Rank Bias Overlap (RBO). So, in the end, we would have a correlation matrix for each similarity matrix where rows and columns are disease- and drug-control samples (or the other way around).

*Note*:
   - Pearson and Spearman assign an equal weight to all genes.
   - RBO allows us to assign higher weights to those are highy expressed (ranked at the top of the list) and lower weights to those with low expression values i.e., noises.
   - Make sure to have all the genes are in the same order when computing correlation.
   - In our analysis, there is no a specific cutoff for correlation to say that a correlation value is meaningful since we are working with quite big lists with ~1000 genes meaning that it's unlikely to have a very high correlation close to one (or it could be..). In practice, a correlation of 0.7 or 0.5 could probably be meaningful too if they happen to be higher than the resulting average correlation.
5. Result interpretation
   - Normalize the correlation matrices where both row and column vectors are taken into account (mathematical explaination will be discussed later), so that we can compare the correlation across all pairs.
     - Use z-score normalization $z_{i,j} = (z_i + z_j)/\sqrt{2}$ where $z_i = (r_{i,j} - \mu_{i})/\sigma_{i}$, $z_j = (r_{i,j} - \mu_{j})/\sigma_{j}$, $r_{i,j}$ is the correlation value at row $i$ column $j$.
   - We should be able to see some clusters of sample pairs with high correlation that are originally coming from (biologically/morphologically) same/similar tissues and cell lines.
   - For each disease control sample (on the row), sort the drug control samples by their correlation values (column sorting).
6. Result validation
   - Rerun the sample analysis, but this time shuffle the expression values so that we have samples with the sample genes but with random expression values. We do this randomization to check if the results we got from the analysis using real data is meaningful; if yes, we should not see any significant pattern in our random-expression results i.e., the correlation should be close to zero.
