# TB data currations documentations

Author : Ling Thang
Date : 2026-03-13

Preprocessed transcriptomic datasets were retrieved from ARCHS4 and their respective metadata from Gene Expression Omnibus (GEO) databases, followed by rigorous data curation, preprocessing, and differential expression analysis to identify disease signatures.

## Data curation and preprocessing

### [Stage1.0_getDiseases-ling.ipynb](Stage1.0_getDiseases-ling.ipynb)

Querying ARCHS4 for disease-specific transcriptomic datasets and retrieving metadata from GEO.

- for full reproduciblity please install the files from archs4 listed in the notebook

**Steps after this notebook include :**

- querying Gene Expression Omnibus (GEO) for the metadata,
  - see [get_metadata_GEO.R](get_metadata_GEO.R) and [get_expression_data_ARCHS4.R](get_expression_data_ARCHS4.R) for retrieval script
- manually identifying relevant samples for for the differential expression analysis,
  - see [data/labeling/TB_currationdoc.tsv](https://github.com/JRaviLab/integrative-drugrep-tb/blob/main/data/labeling/TB_currationdoc.tsv) for specifc details of justification for inclusion/exclusion of samples

### [Stage1.1_CheckLabeling.ipynb](Stage1.1_CheckLabeling.ipynb)

Verifying the accuracy of manual labeling of samples based on metadata and curated lists.

- standardizes the cell and tissues sources in the geo metadata

### [Stage1.2_cleanup_expression_data.ipynb](Stage1.2_cleanup_expression_data.ipynb)

- please unzip the microarray raw expression files 1 and 2 and place them in a folder named `microarray_data_forDE/rawexpression`
- Ensures that of the samples in the GEO metadata are present in the Archs4 expression data

### [Stage1.3_CriteriaCheck.ipynb](Stage1.3_CriteriaCheck.ipynb)

Checking metadata criteria for DESEQ.

- Ensures that every study has at least 3 samples in both disease and control groups
