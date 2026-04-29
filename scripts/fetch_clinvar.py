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


def esearch_clinvar(gene: str, api_key: str) -> list[str]:
    """Search ClinVar for missense variants of a gene. Return list of UIDs."""
    term = f'"{gene}"[Gene] AND single_gene[prop] AND "missense variant"[molecular_consequence]'
    params = {
        "db": "clinvar",
        "term": term,
        "retmax": 5000,
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key

    url = f"{EUTILS_BASE}/esearch.fcgi"
    log.info(f"esearch: {url} term={term}")
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    count = int(data["esearchresult"]["count"])
    ids = data["esearchresult"]["idlist"]
    log.info(f"Found {count} ClinVar records, retrieved {len(ids)} IDs")
    return ids


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
        resp = requests.get(url, params=params)
        resp.raise_for_status()
        data = resp.json()

        for uid in batch:
            if uid in data.get("result", {}):
                all_records.append(data["result"][uid])

        time.sleep(0.35 if not api_key else 0.11)
    return all_records


def parse_cdna_position(cdna_change: str) -> int | None:
    """Extract mRNA position from cDNA notation like 'c.761G>A' -> 761."""
    if not cdna_change:
        return None
    m = re.match(r"c\.(\d+)", cdna_change)
    return int(m.group(1)) if m else None


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

    mrna_pos = parse_cdna_position(cdna_change)

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
        "mrna_position": mrna_pos,
    }


def main():
    gene = snakemake.params.gene
    api_key = snakemake.params.api_key
    output_mutations = snakemake.output.mutations
    output_transcript = snakemake.output.transcript_id

    log.info(f"Fetching ClinVar missense mutations for gene: {gene}")

    ids = esearch_clinvar(gene, api_key)
    if not ids:
        log.error(f"No ClinVar records found for {gene}")
        sys.exit(1)

    records = esummary_batch(ids, api_key)
    log.info(f"Fetched {len(records)} summary records")

    parsed = []
    transcript_counter = Counter()
    for rec in records:
        p = parse_record(rec)
        if p and p["mrna_position"] is not None:
            parsed.append(p)
            transcript_counter[p["transcript_id"]] += 1

    log.info(f"Parsed {len(parsed)} records with valid mRNA positions")

    if not parsed:
        log.error("No valid missense mutation records found")
        sys.exit(1)

    principal_transcript = transcript_counter.most_common(1)[0][0]
    log.info(
        f"Principal transcript: {principal_transcript} "
        f"({transcript_counter[principal_transcript]}/{len(parsed)} records)"
    )

    df = pd.DataFrame(parsed)
    df = df.sort_values("mrna_position").reset_index(drop=True)
    df.to_csv(output_mutations, sep="\t", index=False)
    log.info(f"Wrote {len(df)} mutations to {output_mutations}")

    Path(output_transcript).write_text(principal_transcript + "\n")
    log.info(f"Wrote principal transcript to {output_transcript}")


if __name__ == "__main__":
    main()
