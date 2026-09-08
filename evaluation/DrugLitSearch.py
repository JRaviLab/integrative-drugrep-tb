# Literature searcher for BROAD's list of FDA approved drugs
# original author : @LingT03
# last modified : 09/03/2026
# LLM usage note : 
# The template was generated with the help of an LLM and then reviewed 
# and modified to fit the project's needs by the author. Systematic code reviews were 
# performed by manuscript co-authors to ensure the accuracy of the script and it's outputs
#
# Each retrieved record is assigned one evidence level from MEDLINE
# publication types and MeSH headings:
#
#   clinical_study > human_subject > animal > in_vitro > other
#
# The levels are mutually exclusive, so the per-drug counts sum to the
# number of primary records for that drug. Reviews are counted separately
# and excluded from the levels

import argparse
import configparser
import json
import os
import re
from collections import Counter
from typing import Iterable, Iterator

import pandas as pd
from Bio import Entrez, Medline
from tqdm import tqdm

# config
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_config = configparser.ConfigParser()
_config.read(os.path.join(_SCRIPT_DIR, '.config'))

Entrez.email = _config.get('DEFAULT', 'email')
Entrez.api_key = _config.get('DEFAULT', 'api_key', fallback=None)

# parameters
MAX_RESULTS_PER_DRUG = 1000
EFETCH_BATCH_SIZE    = 200
CACHE_VERSION        = 2

# queries

# MeSH-anchored TB terms (papers about TB), explosion left on so
# extrapulmonary forms are included
TB_MESH = ('"Tuberculosis"[MeSH Terms] OR '
           '"Tuberculosis, Pulmonary"[MeSH Terms] OR '
           '"Mycobacterium tuberculosis"[MeSH Terms]')

# free-text fallback for papers not yet MeSH-indexed
TB_FREE = ('"Mycobacterium tuberculosis"[Title/Abstract] OR '
           '"tuberculosis"[Title/Abstract]')

DISEASE_TERMS = f'({TB_MESH} OR {TB_FREE})'

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


def build_drug_query(drug_name: str) -> str:
    """
    One query per drug, drug text or exact MeSH descriptor, anchored to TB

    :noexp keeps a descriptor from exploding onto its children
    """
    safe = drug_name.replace('"', "")
    drug_terms = f'("{safe}"[Title/Abstract] OR "{safe}"[MeSH Terms:noexp])'
    return f"({drug_terms}) AND ({DISEASE_TERMS})"


# evidence definitions

# PT is a closed NLM vocabulary
CLINICAL_PT = frozenset({
    "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial", "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial", "Equivalence Trial",
    "Observational Study", "Multicenter Study",
})
CASE_PT   = frozenset({"Case Reports"})
REVIEW_PT = frozenset({"Review", "Systematic Review", "Meta-Analysis"})

# the Humans check tag also covers cultured human cells, so subject-level
HUMAN_SUBJECT_MH = frozenset({
    "Adult", "Middle Aged", "Aged", "Aged, 80 and over", "Young Adult",
    "Adolescent", "Child", "Child, Preschool", "Infant", "Infant, Newborn",
    "Cohort Studies", "Prospective Studies", "Retrospective Studies",
    "Case-Control Studies", "Cross-Sectional Studies", "Treatment Outcome",
})

ANIMAL_MH = frozenset({
    "Animals", "Disease Models, Animal", "Mice", "Rats", "Guinea Pigs",
    "Rabbits", "Zebrafish", "Macaca", "Macaca mulatta", "Macaca fascicularis",
})

# culture context headings only
IN_VITRO_MH = frozenset({
    "In Vitro Techniques", "Cell Culture Techniques", "Tissue Culture Techniques",
    "Cell Culture Techniques, Three Dimensional", "Primary Cell Culture",
    "Coculture Techniques", "Organoids",
    "Cells, Cultured", "Cell Line", "Cell Line, Tumor",
    "THP-1 Cells", "A549 Cells",
})

# fallback keywords for records with no MeSH indexing
HUMAN_SUBJECT_KW = ("patients", "participants", "clinical trial", "randomized",
                    "cohort study", "case-control")
ANIMAL_KW        = ("mice", "mouse model", "murine", "guinea pig", "macaque",
                    "zebrafish", "rabbit", "rat model")
IN_VITRO_KW      = ("in vitro", "cell culture", "cultured cells", "cell line",
                    "organoid", "thp-1")

# ordered strongest to weakest
# 'other' holds work unrelated to HDT evidence
EVIDENCE_LEVELS = ("clinical_study", "human_subject", "animal", "in_vitro", "other")

# classification

def _as_text(value: object) -> str:
    """Flatten a MEDLINE field to a string, OAB and OT are list-valued"""
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)):
        return " ".join(_as_text(v) for v in value)
    return "" if value is None else str(value)


def _normalize_mesh(mesh_terms: Iterable[str]) -> set[str]:
    """Strip major-topic asterisks and subheadings from raw MH entries"""
    return {t.lstrip("*").split("/")[0].strip() for t in mesh_terms if t}


def _record_text(record: dict) -> str:
    """Everything [Title/Abstract] covers, tiab includes OT and OAB"""
    parts = (record.get("title", ""), record.get("abstract", ""),
             record.get("other_abstract", ""), record.get("keywords", []))
    return " ".join(t for t in (_as_text(p) for p in parts) if t)


def _matches(text: str, terms: Iterable[str]) -> bool:
    """Word-boundary phrase match, keeps 'rat' out of 'strategy'"""
    return any(re.search(rf"(?<!\w){re.escape(t)}(?!\w)", text, re.I) for t in terms)


def _level_from_mesh(pts: set[str], mesh: set[str]) -> str:
    """Assign an evidence level from publication types and MeSH headings"""
    in_vitro = bool(mesh & IN_VITRO_MH)

    if pts & CLINICAL_PT:
        return "clinical_study"
    if "Humans" in mesh and (mesh & HUMAN_SUBJECT_MH or pts & CASE_PT):
        return "human_subject"
    # the Animals tag alone can be an animal-derived cell line
    if mesh & ANIMAL_MH and ("Disease Models, Animal" in mesh or not in_vitro):
        return "animal"
    if in_vitro:
        return "in_vitro"
    return "other"


def _level_from_text(record: dict) -> str:
    """Lower-confidence fallback for records with no MeSH indexing"""
    text = _record_text(record)
    if _matches(text, HUMAN_SUBJECT_KW):
        return "human_subject"
    if _matches(text, ANIMAL_KW) and not _matches(text, IN_VITRO_KW):
        return "animal"
    if _matches(text, IN_VITRO_KW):
        return "in_vitro"
    return "other"


def classify_record(record: dict) -> dict:
    """
    Assign one evidence level and a review flag to a single record

    Records without MeSH indexing (preprints, PMC-only) fall back to
    title/abstract/keyword matching and carry mesh_indexed=False so the
    lower-confidence subset stays auditable
    """
    pts     = set(record.get("pub_types", []))
    mesh    = _normalize_mesh(record.get("mesh", []))
    indexed = bool(mesh)

    return {
        "evidence_level": _level_from_mesh(pts, mesh) if indexed
                          else _level_from_text(record),
        # a case report that also reviews the literature is still primary
        "is_review": bool(pts & REVIEW_PT) and not bool(pts & (CLINICAL_PT | CASE_PT)),
        "mesh_indexed": indexed,
    }


def summarize_drug(records: list[dict]) -> dict:
    """
    Per-drug evidence counts

    Reviews are counted separately and excluded from the levels, so the
    n_* level counts sum to the number of primary records
    """
    primary = [r for r in records if not r["is_review"]]
    levels  = Counter(r["evidence_level"] for r in primary)

    counts: dict[str, object] = {f"n_{lvl}": levels[lvl] for lvl in EVIDENCE_LEVELS}
    counts["n_review"]      = sum(r["is_review"] for r in records)
    counts["n_not_indexed"] = sum(not r["mesh_indexed"] for r in primary)
    counts["n_hdt"]         = sum(r["is_hdt"] for r in records)
    counts["highest_evidence_level"] = next(
        (lvl for lvl in EVIDENCE_LEVELS if levels[lvl]), "none"
    )
    return counts


# searching

def canonical_key(pmid: str = "", pmc_id: str = "") -> str:
    """
    Stable per-article key

    a PMID and a PMCID name the same article with different strings, so
    counting or deduplicating on the raw ids double counts PubMed/PMC twins
    PMID wins because PubMed is the citation source
    """
    pmid, pmc_id = str(pmid).strip(), str(pmc_id).strip()
    if pmid:
        return f"pmid:{pmid}"
    return f"pmcid:{pmc_id}" if pmc_id else ""


def _batched(items: list, size: int) -> Iterator[list]:
    """Yield fixed size chunks, keeps efetch requests manageable"""
    for start in range(0, len(items), size):
        yield items[start:start + size]


class LiteratureSearcher:
    """
    PubMed (and optionally PMC) literature searcher for TB drug candidates

    Two searches per drug: the unified drug + TB query, and the same query
    restricted to HDT terms. The second exists so HDT papers cannot be
    pushed out of the result set by the retmax cap, and its id set is the
    HDT annotation, PubMed does the matching rather than local text rules

    Bio.Entrez paces requests itself (0.1s with an api key, 0.37s without)
    so there are no manual sleeps here
    """

    def __init__(self, email: str, api_key: str | None = None,
                 search_pmc: bool = False,
                 max_results: int = MAX_RESULTS_PER_DRUG) -> None:
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.search_pmc  = search_pmc
        self.max_results = max_results

    
    def _search(self, db: str, query: str) -> tuple[int, list[str]]:
        try:
            handle = Entrez.esearch(db=db, term=query, retmax=self.max_results)
            record = Entrez.read(handle)
            handle.close()
            return int(record.get("Count", 0)), list(record.get("IdList", []))
        except Exception as e:
            print(f"  - esearch failed on {db}: {e}")
            return 0, []

    def _search_both(self, query: str) -> tuple[int, int, list[str]]:
        """
        Search PubMed and, when enabled, PMC

        Returns (pubmed_count, pmc_count, ids) with the counts kept apart
        because the same article can be hit in both databases, so their sum
        is a hit total and not a study count
        """
        count, ids = self._search("pubmed", query)
        if not self.search_pmc:
            return count, 0, ids

        pmc_count, pmc_ids = self._search("pmc", query)
        prefixed = [i if i.startswith("PMC") else f"PMC{i}" for i in pmc_ids]
        return count, pmc_count, list(dict.fromkeys(ids + prefixed))

    def find_evidence(self, drug_name: str) -> tuple[dict, list[dict]]:
        """
        Returns (counts, summaries) for *drug_name*

        counts holds pubmed_hits, pmc_hits, unique_articles and truncated
        """
        print(f"\nSearching for '{drug_name}'...")

        base = build_drug_query(drug_name)
        pm_count, pmc_count, ids = self._search_both(base)
        truncated = pm_count > self.max_results or pmc_count > self.max_results

        _, _, hdt_ids = self._search_both(f"({base}) AND ({HDT_TERMS})")

        merged = list(dict.fromkeys(ids + hdt_ids))
        counts = {"pubmed_hits": pm_count, "pmc_hits": pmc_count,
                  "unique_articles": 0, "results_truncated": truncated}
        if not merged:
            print("    → no results")
            return counts, []

        cap = f", capped at {self.max_results}" if truncated else ""
        hits = f"{pm_count} PubMed" + (f" + {pmc_count} PMC" if self.search_pmc else "")
        print(f"    → {hits} hits{cap}. Fetching...")

        summaries, pmid_of_pmcid = self._fetch(merged)

        # a PMC hit and its PubMed twin are one article, so resolve every id
        # to a canonical key before testing HDT membership
        hdt_keys = set()
        for i in hdt_ids:
            if str(i).startswith("PMC"):
                hdt_keys.add(canonical_key(pmid_of_pmcid.get(i, ""), i))
            else:
                hdt_keys.add(canonical_key(i))
        for s in summaries:
            s["is_hdt"] = canonical_key(s["pmid"], s["pmc_id"]) in hdt_keys

        counts["unique_articles"] = len(summaries)
        print(f"    → {len(summaries)} unique articles "
              f"({sum(s['is_hdt'] for s in summaries)} HDT)")
        return counts, summaries

    @staticmethod
    def _extract_doi(record: dict) -> str:
        """
        Take the token before a [doi] tag in LID or AID
        LID often holds a [pii] first, so the first token is not safe
        """
        for field in [record.get("LID", "")] + list(record.get("AID", [])):
            tokens = field.split()
            for i, token in enumerate(tokens):
                if token == "[doi]" and i:
                    return tokens[i - 1]
        return ""

    def _fetch(self, ids: list[str]) -> tuple[list[dict], dict[str, str]]:
        """
        Returns (summaries, pmcid -> pmid) with PubMed/PMC twins collapsed
        onto the PubMed record
        """
        pmids   = [i for i in ids if not str(i).startswith("PMC")]
        pmc_ids = [i for i in ids if str(i).startswith("PMC")]

        summaries = self._fetch_medline(pmids) if pmids else []
        if not (pmc_ids and self.search_pmc):
            return summaries, {}

        seen = {s["pmid"] for s in summaries if s["pmid"]}
        pmc_summaries, pmid_of_pmcid = self._fetch_pmc(pmc_ids, seen)

        # record the PMCID on the twin that was already fetched from PubMed
        by_pmid = {s["pmid"]: s for s in summaries if s["pmid"]}
        for pmcid, pmid in pmid_of_pmcid.items():
            if pmid in by_pmid and not by_pmid[pmid]["pmc_id"]:
                by_pmid[pmid]["pmc_id"] = pmcid
        return summaries + pmc_summaries, pmid_of_pmcid

    def _fetch_medline(self, pmids: list[str]) -> list[dict]:
        summaries = []
        for batch in _batched(pmids, EFETCH_BATCH_SIZE):
            try:
                handle = Entrez.efetch(db="pubmed", id=",".join(batch),
                                       rettype="medline", retmode="text")
                for record in Medline.parse(handle):
                    pmid = record.get("PMID", "")
                    if " " in pmid:
                        # a missing separator merged consecutive records
                        print(f"  - warning: merged record skipped ({pmid[:30]}...)")
                        continue
                    summaries.append({
                        "pmid":    pmid,
                        "pmc_id":  "",
                        "doi":     self._extract_doi(record),
                        "title":   record.get("TI", ""),
                        "authors": ", ".join(record.get("AU", [])),
                        "source":  "PubMed",
                        # kept for classification, tiab also covers OT and OAB
                        "pub_types":      list(record.get("PT", [])),
                        "mesh":           list(record.get("MH", [])),
                        "abstract":       record.get("AB", ""),
                        "other_abstract": _as_text(record.get("OAB", "")),
                        "keywords":       list(record.get("OT", [])),
                    })
                handle.close()
            except Exception as e:
                print(f"  - efetch failed for batch of {len(batch)}: {e}")
        return summaries

    def _fetch_pmc(self, pmc_ids: list[str],
                   seen_pmids: set[str]) -> tuple[list[dict], dict[str, str]]:
        """
        Resolve PMC hits via esummary, then route PubMed-linked records
        through the MEDLINE fetch so MeSH and PT are available

        efetch db=pmc rettype=medline is avoided, its output can arrive
        without record separators and Medline.parse then silently merges
        consecutive records into one
        """
        raw_ids  = [i.replace("PMC", "") for i in pmc_ids]
        docsums  = []
        for chunk in _batched(raw_ids, EFETCH_BATCH_SIZE):
            try:
                handle = Entrez.esummary(db="pmc", id=",".join(chunk))
                docsums += Entrez.read(handle)
                handle.close()
            except Exception as e:
                print(f"  - PMC esummary failed: {e}")

        # pmid_of_pmcid covers every resolvable hit, including twins already
        # fetched from PubMed, so HDT keys can be canonicalized later
        linked, pmc_of_pmid, pmid_of_pmcid, unlinked = [], {}, {}, []
        for docsum in docsums:
            ids = docsum.get("ArticleIds", {}) or {}
            if isinstance(ids, list):
                ids = {d.get("idtype", ""): d.get("value", "") for d in ids
                       if isinstance(d, dict)}
            pmcid = str(ids.get("pmcid") or docsum.get("Id", ""))
            if pmcid and not pmcid.startswith("PMC"):
                pmcid = f"PMC{pmcid}"
            pmid = str(ids.get("pmid") or "").strip()

            if pmid and pmid != "0":
                pmid_of_pmcid[pmcid] = pmid
                if pmid not in seen_pmids:
                    linked.append(pmid)
                    pmc_of_pmid[pmid] = pmcid
            else:
                # PMC-only record, esummary metadata is all there is
                unlinked.append({
                    "pmid":    "",
                    "pmc_id":  pmcid,
                    "doi":     str(ids.get("doi") or ""),
                    "title":   str(docsum.get("Title", "")),
                    "authors": ", ".join(str(a) for a in docsum.get("AuthorList", [])),
                    "source":  "PMC",
                    "pub_types": [], "mesh": [], "abstract": "",
                    "other_abstract": "", "keywords": [],
                })

        summaries = self._fetch_medline(list(dict.fromkeys(linked)))
        for s in summaries:
            s["pmc_id"] = pmc_of_pmid.get(s["pmid"], "")
            s["source"] = "PMC"
        return summaries + unlinked, pmid_of_pmcid


# query cache

# jsonl, one completed drug query appended per line, so an interrupted run
# resumes and a classification change can be re-applied without re-querying
def load_cache(path: str) -> dict[str, dict]:
    if not os.path.exists(path):
        return {}
    cache, skipped = {}, 0
    with open(path) as fh:
        for line in fh:
            try:
                entry = json.loads(line)
                if entry["version"] == CACHE_VERSION:
                    cache[entry["key"]] = entry
                else:
                    skipped += 1
            except (json.JSONDecodeError, KeyError):
                skipped += 1
    if skipped:
        print(f"  - ignored {skipped} stale or malformed cache lines")
    return cache


def append_cache(path: str, key: str, counts: dict,
                 summaries: list[dict]) -> None:
    entry = {"version": CACHE_VERSION, "key": key, "counts": counts,
             "summaries": summaries}
    with open(path, "a") as fh:
        fh.write(json.dumps(entry) + "\n")


# output

def format_evidence(summaries: list[dict]) -> str:
    """Human readable evidence block for the TSV"""
    if not summaries:
        return "No evidence found"
    entries = []
    for i, s in enumerate(summaries, 1):
        pmc = f" | PMCID: {s['pmc_id']}" if s["pmc_id"] else ""
        hdt = " | HDT" if s.get("is_hdt") else ""
        entries.append(
            f"[{i}] {s['title']}\n"
            f"    Authors: {s['authors']}\n"
            f"    PMID: {s['pmid']}{pmc} | DOI: {s['doi']} | "
            f"Evidence: {s['evidence_level']}{hdt}"
        )
    return "\n\n".join(entries)

def print_summary(records: dict[str, dict]) -> None:
    """Corpus level distribution, the numbers to quote in Methods"""
    if not records:
        return
    total   = len(records)
    levels  = Counter(r["evidence_level"] for r in records.values())
    reviews = sum(r["is_review"] for r in records.values())
    unindexed = sum(not r["mesh_indexed"] for r in records.values())

    print(f"\nClassified {total} unique records")
    for lvl in EVIDENCE_LEVELS:
        print(f"  {levels[lvl]:6d}  {lvl}")
    print(f"  {reviews:6d}  reviews (counted separately)")
    print(f"  {unindexed} records ({unindexed / total:.1%}) had no MeSH indexing "
          f"and were classified from text")


# CLI

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="DrugLitSearch",
        description="Search PubMed (and optionally PMC) for TB drug literature "
                    "and stratify each record by evidence level.",
    )
    p.add_argument("input", metavar="INPUT_FILE", help="Input CSV or TSV file.")
    p.add_argument("column", metavar="DRUG_COLUMN",
                   help="Column containing drug names.")
    p.add_argument("--article", action="store_true",
                   help="Also search PubMed Central.")
    p.add_argument("--out-tsv", metavar="PATH", default=None,
                   help="Output TSV (default: <input_stem>_literature.tsv).")
    p.add_argument("--out-json", metavar="PATH", default=None,
                   help="Classified records JSON (default: <input_stem>_classified.json).")
    p.add_argument("--cache", metavar="PATH", default=None,
                   help="Query cache (default: <input_stem>_query_cache.jsonl).")
    p.add_argument("--no-cache", action="store_true",
                   help="Ignore and do not write the query cache.")
    p.add_argument("--max-results", metavar="N", type=int,
                   default=MAX_RESULTS_PER_DRUG,
                   help=f"Records per drug (default: {MAX_RESULTS_PER_DRUG}).")
    return p


def read_drug_names(df: pd.DataFrame, column: str) -> list[str]:
    """Split pipe separated cells and deduplicate, preserving order"""
    names, seen = [], set()
    for raw in df[column].dropna().unique():
        for part in (p.strip() for p in str(raw).split("|")):
            if part and part not in seen:
                seen.add(part)
                names.append(part)
    return names


def main() -> None:
    args = build_parser().parse_args()

    if not os.path.isfile(args.input):
        raise SystemExit(f"Error: input file not found: {args.input}")

    sep = "\t" if os.path.splitext(args.input)[1].lower() in {".tsv", ".txt"} else ","
    df_input = pd.read_csv(args.input, sep=sep)
    if args.column not in df_input.columns:
        raise SystemExit(f"Error: column '{args.column}' not found. "
                         f"Available: {', '.join(df_input.columns)}")

    drug_names = read_drug_names(df_input, args.column)

    base       = os.path.splitext(os.path.abspath(args.input))[0]
    out_tsv    = args.out_tsv or base + "_literature.tsv"
    out_json   = args.out_json or base + "_classified.json"
    cache_path = None if args.no_cache else (args.cache or base + "_query_cache.jsonl")

    if os.path.abspath(out_tsv) == os.path.abspath(args.input):
        raise SystemExit(f"Error: --out-tsv would overwrite the input: {args.input}")
    for path in (out_tsv, out_json):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    print(f"Drugs      : {len(drug_names)}")
    print(f"Database(s): {'PubMed + PMC' if args.article else 'PubMed'}")
    print(f"NCBI key   : {'present' if Entrez.api_key else 'missing (3 req/sec cap)'}")
    print(f"Output TSV : {out_tsv}")
    print(f"Output JSON: {out_json}")

    cache = load_cache(cache_path) if cache_path else {}
    print(f"Cache      : {cache_path or 'disabled'} ({len(cache)} drugs)\n")

    searcher = LiteratureSearcher(Entrez.email, Entrez.api_key,
                                  args.article, args.max_results)
    rows: list[dict] = []
    records: dict[str, dict] = {}  # record id → classified record, deduped across drugs

    for drug in tqdm(drug_names, desc="Processing drugs"):
        key = drug.lower()
        if key in cache:
            counts, summaries = cache[key]["counts"], cache[key]["summaries"]
        else:
            counts, summaries = searcher.find_evidence(drug)
            if cache_path:
                append_cache(cache_path, key, counts, summaries)

        # classify each record once and reuse it across drugs
        drug_records = []
        for s in summaries:
            rid = canonical_key(s["pmid"], s["pmc_id"])
            if not rid:
                continue
            if rid not in records:
                records[rid] = {**s, **classify_record(s), "drugs": []}
            record = records[rid]
            record["is_hdt"] |= bool(s.get("is_hdt"))
            if drug not in record["drugs"]:
                record["drugs"].append(drug)
            s["evidence_level"] = record["evidence_level"]
            drug_records.append(record)

        rows.append({
            "drug_name": drug,
            # pubmed_hits and pmc_hits are per database and can overlap,
            # unique_articles is the deduplicated denominator to report
            **counts,
            **summarize_drug(drug_records),
            "pmid":                [s["pmid"] for s in summaries],
            "literature_evidence": summaries,
        })

    df = pd.DataFrame(rows)
    df["literature_evidence"] = df["literature_evidence"].apply(format_evidence)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nDone. TSV  → {out_tsv}")

    # the json holds the deduplicated classified records, complementary to the tsv
    with open(out_json, "w") as fh:
        json.dump(list(records.values()), fh, indent=4)
    print(f"       JSON → {out_json} ({len(records)} unique records)")

    print_summary(records)


if __name__ == "__main__":
    main()