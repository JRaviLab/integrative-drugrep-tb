# TB data currations documentations

Author : Ling Thang
Date : 2026-03/13

Preprocessed transcriptomic datasets were retrieved from ARCHS4 and their respective metadata from Gene Expression Omnibus (GEO) databases, followed by rigorous data curation, preprocessing, and differential expression analysis to identify disease signatures.

## Data curation and preprocessing

[Data Curration](S0_getDiseases-ling.ipynb) - Querying ARCHS4 for disease-specific transcriptomic datasets and retrieving metadata from GEO.
[Check Manual Labeling](S1.1_CheckLabeling.ipynb) - Verifying the accuracy of manual labeling of samples based on metadata and curated lists.
[Criteria Check for DESEQ](S1.2_CriteriaCheck.ipynb) - Checking metadata criteria for DESEQ.
[Expression Data Processing](S1.3_expressionDataProcessing.ipynb) - Processing expression data for DESEQ.
[Make Contrast](S1.5_MakeContrast.ipynb) - Creating contrasts for DESEQ.
[DESEQ Analysis](Stage2_VerifyingDGEResults.ipynb) - Performing DESEQ analysis and verifying results for TB.
