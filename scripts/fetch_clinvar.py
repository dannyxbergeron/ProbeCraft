"""Fetch ClinVar missense mutations for a gene via NCBI E-utilities."""

import json
import logging
import re
import sys
import time
from collections import Counter
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 200
MAX_RETRIES = 5


def build_consequence_filter(consequences: str | list[str]) -> str:
    """Build the molecular_consequence portion of a ClinVar esearch term.

    Args:
        consequences: "all" to skip filtering, or a list of consequence strings.

    Returns:
        Empty string if "all", else an AND clause like:
        AND ("type1"[molecular_consequence] OR "type2"[molecular_consequence])
    """
    if consequences == "all":
        return ""

    if isinstance(consequences, str):
        consequences = [consequences]

    if len(consequences) == 1:
        return f' AND "{consequences[0]}"[molecular_consequence]'

    parts = " OR ".join(
        f'"{c}"[molecular_consequence]' for c in consequences
    )
    return f" AND ({parts})"


def _get_delay(api_key: str) -> float:
    """Return inter-request delay. Without API key, use a longer delay
    to stay within NCBI's 3 req/s limit even with parallel jobs."""
    return 0.11 if api_key else 0.5


def _rate_limited_get(url: str, params: dict, api_key: str) -> requests.Response:
    """GET with retry on 429 (rate limit) using exponential backoff."""
    delay = _get_delay(api_key)
    for attempt in range(MAX_RETRIES):
        time.sleep(delay)
        resp = requests.get(url, params=params)
        if resp.status_code != 429:
            return resp
        wait = (2 ** attempt) * delay
        log.warning(f"429 rate limited, retry {attempt + 1}/{MAX_RETRIES} in {wait:.1f}s")
        time.sleep(wait)
    resp.raise_for_status()
    return resp


def esearch_clinvar(gene: str, api_key: str, consequences: str | list[str]) -> list[str]:
    """Search ClinVar for variants of a gene. Return list of UIDs."""
    term = f'"{gene}"[Gene]{build_consequence_filter(consequences)}'
    params = {
        "db": "clinvar",
        "term": term,
        "retmax": 0,
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    url = f"{EUTILS_BASE}/esearch.fcgi"
    log.info(f"esearch: {url} term={term}")
    resp = _rate_limited_get(url, params, api_key)
    resp.raise_for_status()
    data = resp.json()

    count = int(data["esearchresult"]["count"])
    log.info(f"Found {count} total ClinVar records")

    # Paginate to get all IDs
    all_ids = []
    retmax = 5000
    for retstart in range(0, count, retmax):
        params["retmax"] = retmax
        params["retstart"] = retstart
        resp = _rate_limited_get(url, params, api_key)
        resp.raise_for_status()
        data = resp.json()
        batch_ids = data["esearchresult"]["idlist"]
        all_ids.extend(batch_ids)

    log.info(f"Retrieved {len(all_ids)} IDs total")
    return all_ids


def esummary_batch(ids: list[str], api_key: str) -> list[dict]:
    """Fetch esummary records in batches of BATCH_SIZE."""
    all_records = []
    for i in range(0, len(ids), BATCH_SIZE):
        batch = ids[i : i + BATCH_SIZE]
        params = {
            "db": "clinvar",
            "id": ",".join(batch),
            "retmode": "json",
        }
        if api_key:
            params["api_key"] = api_key

        url = f"{EUTILS_BASE}/esummary.fcgi"
        log.info(f"esummary batch {i // BATCH_SIZE + 1}: {len(batch)} IDs")
        resp = _rate_limited_get(url, params, api_key)
        resp.raise_for_status()
        data = resp.json()

        for uid in batch:
            if uid in data.get("result", {}):
                all_records.append(data["result"][uid])

    return all_records


def parse_cdna_position(cdna_change: str) -> tuple[int | None, str | None]:
    """Extract cDNA position and type from HGVS notation.

    Handles both short form (c.761G>A) and full HGVS form
    (NM_007055.4(POLR3A):c.-15C>T).

    Returns:
        (position, cdna_type) where:
        - c.761G>A  -> (761, "CDS")
        - c.-5G>A   -> (5, "5UTR")
        - c.*5G>A   -> (5, "3UTR")
        - unparseable -> (None, None)
    """
    if not cdna_change:
        return None, None
    m = re.search(r"c\.-(\d+)", cdna_change)
    if m:
        return int(m.group(1)), "5UTR"
    m = re.search(r"c\.\*(\d+)", cdna_change)
    if m:
        return int(m.group(1)), "3UTR"
    m = re.search(r"c\.(\d+)", cdna_change)
    if m:
        return int(m.group(1)), "CDS"
    return None, None


TRANSCRIPT_RE = re.compile(r"(NM_\d+\.\d+)")


def parse_record(record: dict) -> dict | None:
    """Extract structured fields from a ClinVar esummary record."""
    title = record.get("title", "")
    m = TRANSCRIPT_RE.search(title)
    if not m:
        return None
    transcript_id = m.group(1)

    gene_symbol = record.get("genes", [{}])[0].get("symbol", "")

    variation_set = record.get("variation_set", [{}])
    var = variation_set[0] if variation_set else {}

    protein_change = record.get("protein_change", "") or ""
    cdna_change = var.get("cdna_change", "") or ""
    # Strip leading transcript(gene): prefix (ClinVar includes it for UTR variants)
    cdna_change = re.sub(r'^[A-Z]{2}_\d+\.\d+\([^)]+\):', '', cdna_change)

    clinical_sig = (
        record.get("germline_classification", {}).get("description", "")
        or ""
    )

    var_loc = var.get("variation_loc", [{}])
    loc = var_loc[0] if var_loc else {}
    chrom = loc.get("chr", "")
    pos_start = loc.get("start", "")
    pos_end = loc.get("stop", "")

    accession = record.get("accession", "")

    cdna_pos, cdna_type = parse_cdna_position(cdna_change)

    return {
        "gene": gene_symbol,
        "transcript_id": transcript_id,
        "protein_change": protein_change,
        "cdna_change": cdna_change,
        "clinical_significance": clinical_sig,
        "chrom": chrom,
        "pos_start": pos_start,
        "pos_end": pos_end,
        "clinvar_accession": accession,
        "cdna_position": cdna_pos,
        "cdna_type": cdna_type,
    }


def main():
    gene = snakemake.params.gene
    api_key = snakemake.params.api_key
    consequences = snakemake.params.molecular_consequences
    region = snakemake.params.region
    output_mutations = snakemake.output.mutations
    output_transcript = snakemake.output.transcript_id

    # Auto-expand consequences to include UTR types when region needs them
    if consequences != "all":
        consequences = list(consequences) if isinstance(consequences, list) else [consequences]
        region_lower = region.lower() if region else ""
        needs_5utr = "5utr" in region_lower or region_lower == "all"
        needs_3utr = "3utr" in region_lower or region_lower == "all"
        if needs_5utr and "5 prime UTR variant" not in consequences:
            consequences.append("5 prime UTR variant")
        if needs_3utr and "3 prime UTR variant" not in consequences:
            consequences.append("3 prime UTR variant")

    consequence_desc = "all types" if consequences == "all" else ", ".join(consequences)
    log.info(f"Fetching ClinVar variants for gene: {gene} (consequences: {consequence_desc})")

    ids = esearch_clinvar(gene, api_key, consequences)
    if not ids:
        log.error(f"No ClinVar records found for {gene}")
        sys.exit(1)

    records = esummary_batch(ids, api_key)
    log.info(f"Fetched {len(records)} summary records")

    parsed = []
    dropped_no_transcript = 0
    transcript_counter = Counter()
    for rec in records:
        p = parse_record(rec)
        if not p:
            dropped_no_transcript += 1
            continue
        parsed.append(p)
        transcript_counter[p["transcript_id"]] += 1

    n_with_cdna = sum(1 for p in parsed if p["cdna_position"] is not None)
    log.info(
        f"Parsed {len(parsed)} records ({n_with_cdna} with valid cDNA position, "
        f"{len(parsed) - n_with_cdna} without). "
        f"Dropped {dropped_no_transcript} records (no NM_ transcript in title)."
    )

    if not parsed:
        log.error("No valid missense mutation records found")
        sys.exit(1)

    principal_transcript = transcript_counter.most_common(1)[0][0]
    log.info(
        f"Principal transcript: {principal_transcript} "
        f"({transcript_counter[principal_transcript]}/{len(parsed)} records)"
    )

    # Sort by cdna_type then cdna_position (records without go to the end)
    type_order = {"5UTR": 0, "CDS": 1, "3UTR": 2}
    df = pd.DataFrame(parsed)
    df["_sort_key"] = df["cdna_type"].map(type_order).fillna(3)
    df = df.sort_values(
        ["_sort_key", "cdna_position"], na_position="last"
    ).reset_index(drop=True)
    df = df.drop(columns=["_sort_key"])
    df.to_csv(output_mutations, sep="\t", index=False)
    log.info(f"Wrote {len(df)} mutations to {output_mutations}")

    Path(output_transcript).write_text(principal_transcript + "\n")
    log.info(f"Wrote principal transcript to {output_transcript}")


if __name__ == "__main__":
    main()
