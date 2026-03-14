# TB data currations documentations

Author : Ling Thang
Date : 2026-03-13

Preprocessed transcriptomic datasets were retrieved from ARCHS4 and their respective metadata from Gene Expression Omnibus (GEO) databases, followed by rigorous data curation, preprocessing, and differential expression analysis to identify disease signatures.

## Data curation and preprocessing

[Get Disease](S1.0_getDiseases-ling.ipynb) - Querying ARCHS4 for disease-specific transcriptomic datasets and retrieving metadata from GEO.

- for full reproduciblity please install the files from archs4 listed in the notebook
- this was not included in the data folder as the combined size of the two files are over 50 GB

[Check Manual Labeling](S1.1_CheckLabeling.ipynb) - Verifying the accuracy of manual labeling of samples based on metadata and curated lists.

- standardizes the cell and tissues sources in the geo metadata

[Expression Data Clean Up](Stage1.2_cleanup_expression_data.ipynb)

- Ensures that of the samples in the GEO metadata are present in the Archs4 expression data

[Checking Metadata Criteria](S1.3_CriteriaCheck.ipynb) - Checking metadata criteria for DESEQ.

- Ensures that every study has at least 3 samples in both disease and control groups

[Make Contrasts for DESEQ](Stage1.4_MakeContrast.ipynb)

- This contains the logic for creating the contrast for DESEQ as well as the input file used for the Differential Gene Expression (DGE) analysis.

- contrast are a unique combination of the following metadata fields : series_id, UBERON:NAME, CL:NAME, CONTEXT

- additionally ensure thats all contrasts have at least 3 samples in both disease and control groups
