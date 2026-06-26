import pandas as pd
import requests
import time
import argparse
import os
from Bio import Entrez, Medline
from tqdm import tqdm
import configparser

# config
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_config = configparser.ConfigParser()
_config.read(os.path.join(_SCRIPT_DIR, '.config'))

Entrez.email = _config.get('DEFAULT', 'email')
Entrez.api_key = _config.get('DEFAULT', 'api_key', fallback=None)

# parameters
MAX_RESULTS_PER_DRUG = 1000

# query term constants
# MeSH-anchored TB terms (papers *about* TB, not just mentioning it)
TB_MESH = ('"Tuberculosis"[MeSH Terms] OR '
           '"Tuberculosis, Pulmonary"[MeSH Terms] OR '
           '"Mycobacterium tuberculosis"[MeSH Terms]')

# Free-text fallback for papers not yet MeSH-indexed
TB_FREE = ('"Mycobacterium tuberculosis"[Title/Abstract] OR '
           '"tuberculosis"[Title/Abstract]')

DISEASE_TERMS = f'({TB_MESH} OR {TB_FREE})'

# HDT terms anchored to Title/Abstract to avoid matching footnotes/acknowledgements
HDT_TERMS = (
    '"host-directed therapy"[Title/Abstract] OR '
    '"host-directed therapies"[Title/Abstract] OR '
    '"host-directed drug"[Title/Abstract] OR '
    '"host-directed drugs"[Title/Abstract] OR '
    '"adjunctive therapy"[Title/Abstract] OR '
    '"adjunctive therapies"[Title/Abstract] OR '
    '"adjunctive therapeutic"[Title/Abstract] OR '
    '"adjunctive therapeutics"[Title/Abstract] OR '
    '"host-targeting therapy"[Title/Abstract] OR '
    '"host-targeting therapies"[Title/Abstract] OR '
    '"HDT"[Title/Abstract]'
)

PA_TAG = "[Pharmacological Action]"


# lit searcher class
class LiteratureSearcher:
    """
    Tiered PubMed (and optionally PMC) literature searcher for TB drug candidates.

    Default search tiers per drug:
      1.  Drug name text  + TB + HDT        (direct, specific)
      2.  Drug name text  + TB              (direct, broad)
      3.  Drug name MeSH  + TB + HDT        (MeSH descriptor, specific)
      4.  Drug name MeSH  + TB              (MeSH descriptor, broad)
      5.  Gene targets via DGIdb + TB + HDT (mechanistic, specific)
      6.  Gene targets via DGIdb + TB       (mechanistic, broad)

    With --pa:
      PA tiers are inserted between MeSH (4) and DGIdb (5/6):
      5.  Drug name PA tag + TB + HDT
      6.  Drug name PA tag + TB
      7/8. DGIdb (as above, renumbered)

    With --article: each tier also queries PMC; results are deduplicated.
    """

    def __init__(self, email, api_key=None, search_pmc=False, use_pa=False):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.search_pmc = search_pmc
        self.use_pa     = use_pa
        self.dgidb_url  = "https://dgidb.org/api/v2/interactions.json"

    # Field-scope heuristic
    def _is_specific_drug_name(self, drug_name: str) -> bool:
        """
        Short, specific names (e.g. 'metformin') → [All Fields].
        Long mechanistic class names (e.g. 'prostanoid receptor antagonist')
        → [Title/Abstract] to avoid noise.
        """
        if not isinstance(drug_name, str):
            return False
        generic_kw = {
            'receptor', 'inhibitor', 'agonist', 'antagonist',
            'enzyme', 'kinase', 'channel', 'transporter', 'synthase',
        }
        words = drug_name.lower().split()
        return len(words) <= 2 and not any(w in generic_kw for w in words)

    # Search helpers
    def _search_db(self, db: str, query: str):
        try:
            handle = Entrez.esearch(db=db, term=query, retmax=MAX_RESULTS_PER_DRUG)
            record = Entrez.read(handle)
            handle.close()
            return int(record.get('Count', 0)), record.get('IdList', [])
        except Exception as e:
            print(f"  - esearch failed on {db} for '{query[:60]}...': {e}")
            return 0, []

    def _search_pubmed(self, query: str):
        return self._search_db("pubmed", query)

    def _search_pmc(self, query: str):
        return self._search_db("pmc", query)

    def _combined_search(self, query: str):
        """
        Search PubMed and (when enabled) PMC. Returns merged deduplicated
        (count, id_list). PMC IDs are prefixed with 'PMC' to stay distinct
        from plain PMIDs.
        """
        count, pmids = self._search_pubmed(query)
        time.sleep(0.34)

        if self.search_pmc:
            pmc_count, pmc_ids = self._search_pmc(query)
            time.sleep(0.34)
            pmc_ids_prefixed = [
                f"PMC{i}" if not i.startswith("PMC") else i for i in pmc_ids
            ]
            merged_ids   = list(dict.fromkeys(pmids + pmc_ids_prefixed))
            merged_count = count + pmc_count
            if pmc_count:
                print(f"    PubMed: {count} | PMC: {pmc_count} → merged: {len(merged_ids)}")
            return merged_count, merged_ids

        return count, pmids

    # Fetch helpers
    def _fetch_pubmed_summaries(self, ids: list) -> list:
        if not ids:
            return []
        plain_pmids = [i for i in ids if not str(i).startswith("PMC")]
        pmc_ids     = [i for i in ids if str(i).startswith("PMC")]
        summaries   = []
        if plain_pmids:
            summaries.extend(self._fetch_medline(plain_pmids))
        if pmc_ids and self.search_pmc:
            summaries.extend(self._fetch_pmc_summaries(pmc_ids))
        return summaries

    def _fetch_medline(self, pmids: list) -> list:
        summaries = []
        try:
            handle  = Entrez.efetch(db="pubmed", id=",".join(pmids),
                                    rettype="medline", retmode="text")
            records = Medline.parse(handle)
            for record in records:
                doi = ""
                if "LID" in record and "[doi]" in record["LID"]:
                    doi = record["LID"].split(" ")[0]
                else:
                    for aid in record.get("AID", []):
                        if "[doi]" in aid:
                            doi = aid.split(" ")[0]
                            break
                summaries.append({
                    "title":   record.get("TI", "No Title Found"),
                    "authors": ", ".join(record.get("AU", [])),
                    "pmid":    record.get("PMID", ""),
                    "pmc_id":  "",
                    "doi":     doi,
                    "source":  "PubMed",
                })
            handle.close()
        except Exception as e:
            print(f"  - PubMed efetch failed: {e}")
        return summaries

    def _fetch_pmc_summaries(self, pmc_ids: list) -> list:
        summaries = []
        raw_ids = [i.replace("PMC", "") for i in pmc_ids]
        try:
            handle  = Entrez.efetch(db="pmc", id=",".join(raw_ids),
                                    rettype="medline", retmode="text")
            records = Medline.parse(handle)
            for record in records:
                doi = ""
                if "LID" in record and "[doi]" in record["LID"]:
                    doi = record["LID"].split(" ")[0]
                else:
                    for aid in record.get("AID", []):
                        if "[doi]" in aid:
                            doi = aid.split(" ")[0]
                            break
                summaries.append({
                    "title":   record.get("TI", "No Title Found"),
                    "authors": ", ".join(record.get("AU", [])),
                    "pmid":    record.get("PMID", ""),
                    "pmc_id":  record.get("PMC", ""),
                    "doi":     doi,
                    "source":  "PMC",
                })
            handle.close()
        except Exception as e:
            print(f"  - PMC efetch failed: {e}")
        return summaries

    def _query_dgidb(self, drug_name: str) -> list:
        params = {'drugs': drug_name, 'interaction_sources': 'DrugBank'}
        try:
            resp = requests.get(self.dgidb_url, params=params, timeout=10)
            resp.raise_for_status()
            data  = resp.json()
            genes = {
                item['geneName']
                for item in data.get('matchedTerms', [])[0].get('interactions', [])
            }
            return list(genes)
        except (requests.exceptions.RequestException, IndexError) as e:
            print(f"  - DGIdb query failed for '{drug_name}': {e}")
            return []

    # Tier runner─
    def _run_tier(self, label: str, query: str):
        """
        Run a single tier. Returns (search_type, count, summaries) if results
        are found, otherwise None.
        """
        print(f"  - {label}...")
        count, ids = self._combined_search(query)
        if ids:
            print(f"    → {count} results. Fetching details...")
            return count, self._fetch_pubmed_summaries(ids)
        return None

    # Main search entry point─
    def find_evidence(self, drug_name: str):
        """
        Run the tiered search for *drug_name*.
        Returns (search_type, count, summaries).
        """
        print(f"\nSearching for '{drug_name}'...")

        field_tag = (
            "[All Fields]" if self._is_specific_drug_name(drug_name)
            else "[Title/Abstract]"
        )
        print(f"  - Field scope: {field_tag} | PMC: {'on' if self.search_pmc else 'off'} | PA: {'on' if self.use_pa else 'off'}")

        # Tiers 1 & 2: direct text search
        result = self._run_tier(
            "Tier 1: text + TB + HDT",
            f'("{drug_name}"{field_tag}) AND ({DISEASE_TERMS}) AND ({HDT_TERMS})',
        )
        if result:
            return "Direct (HDT)", *result

        result = self._run_tier(
            "Tier 2: text + TB (broad)",
            f'("{drug_name}"{field_tag}) AND ({DISEASE_TERMS})',
        )
        if result:
            return "Direct (TB)", *result

        # Tiers 3 & 4: MeSH descriptor search
        result = self._run_tier(
            "Tier 3: MeSH descriptor + TB + HDT",
            f'"{drug_name}"[MeSH Terms] AND ({DISEASE_TERMS}) AND ({HDT_TERMS})',
        )
        if result:
            return "MeSH (HDT)", *result

        result = self._run_tier(
            "Tier 4: MeSH descriptor + TB (broad)",
            f'"{drug_name}"[MeSH Terms] AND ({DISEASE_TERMS})',
        )
        if result:
            return "MeSH (TB)", *result

        # Tiers 5 & 6 (optional): MeSH pharmacological action
        if self.use_pa:
            result = self._run_tier(
                "Tier 5: MeSH PA + TB + HDT",
                f'"{drug_name}"{PA_TAG} AND ({DISEASE_TERMS}) AND ({HDT_TERMS})',
            )
            if result:
                return "MeSH PA (HDT)", *result

            result = self._run_tier(
                "Tier 6: MeSH PA + TB (broad)",
                f'"{drug_name}"{PA_TAG} AND ({DISEASE_TERMS})',
            )
            if result:
                return "MeSH PA (TB)", *result

        # Final tiers: mechanistic via DGIdb
        print("  - Mechanistic: querying DGIdb for gene targets...")
        targets = self._query_dgidb(drug_name)
        time.sleep(0.34)

        if not targets:
            print("    → No gene targets found in DGIdb.")
            return "No Results", 0, []

        print(f"    → Targets: {', '.join(targets[:5])}{'...' if len(targets) > 5 else ''}")
        # Gene names are short and ambiguous — always anchor to Title/Abstract
        target_query_part = " OR ".join(
            [f'"{t}"[Title/Abstract]' for t in targets]
        )

        result = self._run_tier(
            "Mechanistic: gene targets + TB + HDT",
            f'({target_query_part}) AND ({DISEASE_TERMS}) AND ({HDT_TERMS})',
        )
        if result:
            return "Mechanistic (HDT)", *result

        result = self._run_tier(
            "Mechanistic: gene targets + TB (broad)",
            f'({target_query_part}) AND ({DISEASE_TERMS})',
        )
        if result:
            return "Mechanistic (TB)", *result

        print("    → No results found.")
        return "No Results", 0, []


# Output formatter
def format_evidence_for_tsv(summaries: list) -> str:
    if not summaries or not isinstance(summaries, list):
        return "No evidence found"
    entries = []
    for i, s in enumerate(summaries, 1):
        source = f" [{s.get('source', 'PubMed')}]" if s.get('source') else ""
        pmc    = f" | PMCID: {s['pmc_id']}" if s.get('pmc_id') else ""
        entries.append(
            f"[{i}] {s.get('title', 'N/A')}\n"
            f"    Authors: {s.get('authors', 'N/A')}\n"
            f"    PMID: {s.get('pmid', 'N/A')}{pmc} | DOI: {s.get('doi', 'N/A')}{source}"
        )
    return "\n\n".join(entries)


# CLI
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="DrugLitSearch",
        description="Search PubMed (and optionally PMC) for TB host-directed therapy evidence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # PubMed only, moa column
  python DrugLitSearch.py data/launch.csv moa

  # PubMed + PMC, with MeSH PA tiers enabled, custom output paths
  python DrugLitSearch.py data/launch.tsv moa --article --pa \\
      --out-tsv results/lit.tsv --out-json results/lit.json
        """,
    )
    p.add_argument("input",  metavar="INPUT_FILE",
                   help="Path to input CSV or TSV file.")
    p.add_argument("column", metavar="DRUG_COLUMN",
                   help="Column name containing drug/MoA names.")
    p.add_argument("--article", action="store_true", default=False,
                   help="Also search PubMed Central (PMC) full-text archive.")
    p.add_argument("--pa", action="store_true", default=False,
                   help="Enable MeSH Pharmacological Action tiers (may increase noise).")
    p.add_argument("--out-tsv",  metavar="PATH", default=None,
                   help="Output TSV path (default: <input_stem>_literature.tsv).")
    p.add_argument("--out-json", metavar="PATH", default=None,
                   help="Output JSON path (default: <input_stem>_literature.json).")
    return p


def default_output_paths(input_path: str):
    """Derive default TSV/JSON paths from the input file path."""
    base = os.path.splitext(os.path.abspath(input_path))[0]
    return base + "_literature.tsv", base + "_literature.json"


# Entry point
if __name__ == "__main__":
    args = build_parser().parse_args()

    # Load input
    if not os.path.isfile(args.input):
        print(f"Error: input file not found: {args.input}")
        raise SystemExit(1)

    ext = os.path.splitext(args.input)[1].lower()
    sep = "\t" if ext in {".tsv", ".txt"} else ","
    df_input = pd.read_csv(args.input, sep=sep)

    if args.column not in df_input.columns:
        available = ", ".join(df_input.columns.tolist())
        print(f"Error: column '{args.column}' not found.\nAvailable columns: {available}")
        raise SystemExit(1)

    # Parse drug names
    # Each raw cell may contain pipe-separated MoA terms and/or an LLM annotation.
    # Format: "term_a|term_b|LLM Annotation : annotated_term"
    # - Pipe-separated terms → individual queries
    # - LLM Annotation marker → query uses annotated_term, display keeps full string
    raw_names  = df_input[args.column].dropna().unique()
    drug_names = []  # list of (query_name, display_name) tuples

    for raw in raw_names:
        parts = [p.strip() for p in raw.split("|") if p.strip()]
        for part in parts:
            if "LLM Annotation :" in part:
                query_name   = part.split("LLM Annotation :")[-1].strip()
                display_name = part  # full string preserved in output
            else:
                query_name   = part
                display_name = part
            drug_names.append((query_name, display_name))

    # Deduplicate on query_name, preserving insertion order
    seen = set()
    drug_names = [
        (q, d) for q, d in drug_names
        if not (q in seen or seen.add(q))
    ]

    print(f"Loaded {len(raw_names)} raw entries → {len(drug_names)} individual queries after splitting.")

    # Resolve output paths
    default_tsv, default_json = default_output_paths(args.input)
    out_tsv  = args.out_tsv  or default_tsv
    out_json = args.out_json or default_json

    # ensure output directories exist
    for path in (out_tsv, out_json):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    db_label = "PubMed + PMC" if args.article else "PubMed"
    print(f"Database(s): {db_label}")
    print(f"PA tiers   : {'enabled' if args.pa else 'disabled'}")
    print(f"Output TSV : {out_tsv}")
    print(f"Output JSON: {out_json}\n")

    # run search
    searcher     = LiteratureSearcher(
        email      = Entrez.email,
        api_key    = Entrez.api_key,
        search_pmc = args.article,
        use_pa     = args.pa,
    )
    results_list = []

    for query_name, display_name in tqdm(drug_names, desc="Processing drugs"):
        search_type, count, summaries = searcher.find_evidence(query_name)
        results_list.append({
            "drug_name":           display_name,
            "search_type":         search_type,
            "studies_found":       count,
            "pmid":                [s.get("pmid", "") for s in summaries],
            "pmc_id":              [s.get("pmc_id", "") for s in summaries],
            "literature_evidence": summaries,
        })

    # Save outputs
    final_df = pd.DataFrame(results_list)

    df_tsv = final_df.copy()
    df_tsv["literature_evidence"] = df_tsv["literature_evidence"].apply(
        format_evidence_for_tsv
    )
    df_tsv.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nDone. TSV  → {out_tsv}")

    final_df.to_json(out_json, orient="records", indent=4)
    print(f"       JSON → {out_json}")