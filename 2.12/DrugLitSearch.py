import pandas as pd
import requests
import time
from Bio import Entrez, Medline
from tqdm import tqdm
import json
import configparser

# # Configuration
configparser = configparser.ConfigParser()
configparser.read('.config')

Entrez.email = configparser.get('DEFAULT', 'email')
Entrez.api_key = configparser.get('DEFAULT', 'api_key', fallback=None)

# search parameters
MAX_RESULTS_PER_DRUG = 1000

# Base query terms for PubMed
DISEASE_TERMS = '"Mycobacterium tuberculosis" OR "tuberculosis" OR "TB"'
HDT_TERMS = '"host-directed therapy" OR "HDT" OR "adjunctive therapy" OR "adjunctive therapeutic" OR "host-targeting therapies"'

# lit search class
class LiteratureSearcher:
    """
    A class to find literature evidence for drug candidates by searching
    PubMed directly and mechanistically via DGIdb.
    """
    def __init__(self, email, api_key=None):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.dgidb_url = "https://dgidb.org/api/v2/interactions.json"

    def _query_dgidb(self, drug_name):
        """Queries DGIdb to find gene targets for a given drug."""
        params = {'drugs': drug_name, 'interaction_sources': 'DrugBank'}
        try:
            response = requests.get(self.dgidb_url, params=params)
            response.raise_for_status()
            data = response.json()
            # extract unique gene names from the interactions
            genes = {item['geneName'] for item in data.get('matchedTerms', [])[0].get('interactions', [])}
            return list(genes)
        except requests.exceptions.RequestException as e:
            print(f"  - DGIdb query failed for '{drug_name}': {e}")
            return []

    def _search_pubmed(self, query):
        """Performs a search on PubMed and returns the count and PMIDs."""
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=MAX_RESULTS_PER_DRUG)
            record = Entrez.read(handle)
            handle.close()
            count = int(record.get('Count', 0))
            pmids = record.get('IdList', [])
            return count, pmids
        except Exception as e:
            print(f"  - PubMed esearch failed for query '{query[:50]}...': {e}")
            return 0, []

    def _fetch_pubmed_summaries(self, pmids):
        """Fetches structured summaries for a list of PMIDs using Medline parsing."""
        if not pmids:
            return []
        
        summaries = []
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="medline", retmode="text")
            # Use Medline.parse for MEDLINE formatted text records
            records = Medline.parse(handle)
            
            for record in records:
                # extract details safely using .get() with default values
                title = record.get("TI", "No Title Found")
                authors = record.get("AU", [])
                pmid = record.get("PMID", "")
                
                # DOI is often in the 'LID' or 'AID' field for articles from pubmed central
                doi = ""
                # check 'LID' first, as it's a common place for the DOI
                if "LID" in record and "[doi]" in record["LID"]:
                    doi = record["LID"].split(" ")[0]
                else:
                    # Fallback to checking 'AID'
                    for aid in record.get("AID", []):
                        if "[doi]" in aid:
                            doi = aid.split(" ")[0]
                            break
                        
                summaries.append({
                    "title": title,
                    "authors": ", ".join(authors),
                    "pmid": pmid,
                    "doi": doi
                })
            handle.close()
            return summaries
        except Exception as e:
            print(f"  - PubMed efetch failed for PMIDs {pmids}: {e}")
            return []

    def find_evidence(self, drug_name):
        """
        Main search function. Implements a tiered search strategy.
        """
        print(f"\nSearching for '{drug_name}'...")
        
        # very specific query including HDT terms
        print("  - Performing direct search (including HDT terms)...")
        specific_direct_query = f'("{drug_name}"[Title/Abstract]) AND ({DISEASE_TERMS}) AND ({HDT_TERMS})'
        count, pmids = self._search_pubmed(specific_direct_query)
        time.sleep(0.34)

        if pmids:
            print(f"  - Found {count} direct results with HDT terms. Fetching details...")
            summaries = self._fetch_pubmed_summaries(pmids)
            return "Direct (HDT)", count, summaries

        # fall back to a broader direct search, if no match for HDT
        print("  - No direct HDT results found. Performing broader direct search...")
        broad_direct_query = f'("{drug_name}"[Title/Abstract]) AND ({DISEASE_TERMS})'
        count, pmids = self._search_pubmed(broad_direct_query)
        time.sleep(0.34)

        if pmids:
            print(f"  - Found {count} direct results. Fetching details...")
            summaries = self._fetch_pubmed_summaries(pmids)
            return "Direct (TB)", count, summaries

        # mechanistic search (if above fails)
        print("  - No direct results found. Performing mechanistic search via DGIdb...")
        targets = self._query_dgidb(drug_name)
        time.sleep(0.34)

        if not targets:
            print("  - No gene targets found in DGIdb.")
            return "No Results", 0, []
            
        print(f"  - Found targets in DGIdb: {', '.join(targets[:3])}...")
        target_query_part = " OR ".join([f'"{target}"[Title/Abstract]' for target in targets])
        
        specific_mechanistic_query = f'({target_query_part}) AND ({DISEASE_TERMS}) AND ({HDT_TERMS})'
        count, pmids = self._search_pubmed(specific_mechanistic_query)
        time.sleep(0.34)
        
        if pmids:
            print(f"  - Found {count} mechanistic results with HDT terms. Fetching details...")
            summaries = self._fetch_pubmed_summaries(pmids)
            return "Mechanistic (HDT)", count, summaries

        # Fallback to broad mechanistic search
        print("  - No mechanistic HDT results found. Performing broader mechanistic search...")
        broad_mechanistic_query = f'({target_query_part}) AND ({DISEASE_TERMS})'
        count, pmids = self._search_pubmed(broad_mechanistic_query)
        time.sleep(0.34)

        if pmids:
            print(f"  - Found {count} mechanistic results. Fetching details...")
            summaries = self._fetch_pubmed_summaries(pmids)
            return "Mechanistic (TB)", count, summaries

        print("  - No mechanistic results found.")
        return "No Results", 0, []

# output
def format_evidence_for_tsv(summaries):
    """
    Formats the list of summaries into a clean, multi-line string for TSV output.
    """
    if not summaries or not isinstance(summaries, list):
        return "No evidence found"
        
    formatted_strings = []
    for i, s in enumerate(summaries):
        # clean up
        entry = (
            f"[{i+1}] {s.get('title', 'N/A')}\n"
            f"    Authors: {s.get('authors', 'N/A')}\n"
            f"    PMID: {s.get('pmid', 'N/A')} | DOI: {s.get('doi', 'N/A')}"
        )
        formatted_strings.append(entry)
        
    return "\n\n".join(formatted_strings)


if __name__ == "__main__":
    # hard coded for now
    try:
        drug_pred_results = pd.read_csv("drug_list.csv")
    
        drug_names = drug_pred_results['drug_name'].unique()

    except FileNotFoundError:
        print("Error: Input file not found. Please check the path.")
        exit()

    # initialize search and output
    searcher = LiteratureSearcher(email=Entrez.email, api_key=Entrez.api_key)
    # This list will store dictionaries, which is easier to work with
    results_list = []

    for drug in tqdm(drug_names, desc="Processing Drugs"):
        search_type, count, summaries = searcher.find_evidence(drug)        
        
        results_list.append({
            "drug_name": drug,
            "search_type": search_type,
            "studies_found": count,
            "pmid" : [s.get('pmid', '') for s in summaries],  # store PMIDs as a list
            "literature_evidence": summaries # store the structured data
        })

    final_df = pd.DataFrame(results_list)

    # hard coded save paths for now - can be made dynamic later    
    df_for_tsv = final_df.copy()
    df_for_tsv['literature_evidence'] = df_for_tsv['literature_evidence'].apply(format_evidence_for_tsv)
    
    output_path_tsv = "results/drug_predictions_with_literature.tsv"
    df_for_tsv.to_csv(output_path_tsv, sep="\t", index=False)
    print(f"\n✅ Literature search complete. Human-readable output saved to: {output_path_tsv}")

    output_path_json = "results/drug_predictions_with_literature.json"
    final_df.to_json(output_path_json, orient="records", indent=4)
    print(f"✅ Structured JSON output saved to: {output_path_json}")