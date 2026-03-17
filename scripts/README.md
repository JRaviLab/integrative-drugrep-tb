# Scripts

Primary codebase for the integrative TB drug repurposing pipeline. Scripts are numbered to reflect execution order.

| Script | Type | Purpose |
|--------|------|---------|
| `00_multids_microarray_DEwithlimma.R` | R | Differential expression across microarray datasets (limma) |
| `00_multids_RNAseq_DEwithDESeq2.R` | R | Differential expression across RNA-seq datasets (DESeq2) |
| `00_background_genes_PA_functions.R` | R | Background gene sets for pathway analysis (GO, KEGG, LINCS) |
| `00_get_crossSigs_enrichedGO.R` | R | GO enrichment across aggregated signatures |
| `01_signature_aggregation_functions.R` | R | Aggregate DE signatures across datasets; compute membership/Jaccard matrices |
| `01_signature_landmark_prep_functions.R` | R | Filter signatures to L1000 landmark genes |
| `02_drugrep_get_prediction_indiv.R`, `02_drugrep_get_prediction_aggr.R` | R | **Main runners** — query signatures against drug databases to get connectivity scores |
| `02_signatureSearch_connectivity_scores_functions.R` | R | CMap 1.0, LINCS, and correlation-based scoring functions |
| `03_summarize_drugs_methodswise_functions.R` | R | Summarize and compare drug predictions across scoring methods |
| `04_data_generation_RankAggregation_functions.R` | R | Prepare drug lists for rank aggregation |
| `04_partial_druglists_RankAggregation_functions.R` | R | Rank aggregation over partial drug lists |
| `04_BiG_RankAggregation_functions.R` | R | BiG (Bayesian-inspired) rank aggregation |
| `04_RankAggregation_utilities.R` | R | Shared utilities for rank aggregation steps |
| `05_high_confidence_drug_prediction.R` | R | Filter and summarize high-confidence drug candidates |
| `05_build_drugtarget_network_functions.R` | R | Build drug–target interaction networks |
| `suppl_baseline_get_drugdata_lincslevel3.py` | Python | Retrieve LINCS level-3 drug data for baseline comparison |
