## info on last run

> last ran : September 25, 2026
>
> ran by : Ling._.T
>
> command used : 
> python3 evaluation/DrugLitSearch.py evaluation/launched_drugs.tsv drug_name \
>  --article \
>  --refresh \
>  --verbose \
>  --out-tsv evaluation/labels/launched_drug_with_literature.tsv \
>  --out-json evaluation/labels/launched_drug_with_literature.json \
>  2>&1 | tee evaluation/labels/run_$(date +%Y%m%d).log
> 

**see run info below**

### Terminal session

```{zsh}
(.venv) ➜  integrative-drugrep-tb git:(evidence_label) ✗ python3 evaluation/DrugLitSearch.py evaluation/launched_drugs.tsv drug_name \
  --article \
  --refresh \
  --verbose \
  --out-tsv evaluation/labels/launched_drug_with_literature.tsv \
  --out-json evaluation/labels/launched_drug_with_literature.json \
  2>&1 | tee evaluation/labels/run_$(date +%Y%m%d).log
Drugs      : 2427
Database(s): PubMed + PMC
NCBI key   : present
Output TSV : evaluation/labels/launched_drug_with_literature.tsv
Output JSON: evaluation/labels/launched_drug_with_literature.json
Cache      : bypassed (--refresh), every drug re-queried


Done. TSV  → evaluation/labels/launched_drug_with_literature.tsv
       JSON → evaluation/labels/launched_drug_with_literature.json (32392 unique records)

Queried 2427 drugs, reused 0 cached results

Classified 32392 unique records
    1686  clinical_study
   11388  human_subject
    3056  animal
    1201  in_vitro
   12378  other
    2683  reviews (counted separately)
  292 positive records (non review HDT)
  4556 records (14.1%) had no MeSH indexing and were classified from text

```