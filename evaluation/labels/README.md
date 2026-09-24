## info on last run

> last ran : September 11, 2026
>
> ran by : Ling._.T
>
> command used python3 evaluation/DrugLitSearch.py evaluation/launched_drugs.tsv drug_name \
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

Processing drugs: 100%|██████████| 2427/2427 [1:02:11<00:00,  1.54s/it]

Done. TSV  → evaluation/labels/launched_drug_with_literature.tsv
       JSON → evaluation/labels/launched_drug_with_literature.json (32330 unique records)

Queried 2427 drugs, reused 0 cached results

Classified 32368 unique records
    1684  clinical_study
   11383  human_subject
    3056  animal
    1201  in_vitro
   12366  other
    2678  reviews (counted separately)
    377   positive records (non review HDT)
  4538 records (14.0%) had no MeSH indexing and were classified from text
```