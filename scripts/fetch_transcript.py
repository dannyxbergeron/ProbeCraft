"""Fetch transcript sequence and exon-exon junctions from NCBI GenBank."""

import logging
import sys
from pathlib import Path

import requests
from Bio import SeqIO
from io import StringIO

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_genbank(transcript_id: str, api_key: str) -> str:
    """Fetch GenBank record for a RefSeq transcript accession."""
    params = {
        "db": "nucleotide",
        "id": transcript_id,
        "rettype": "gb",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key

    url = f"{EUTILS_BASE}/efetch.fcgi"
    log.info(f"efetch GenBank for {transcript_id}")
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    return resp.text


def parse_cds_feature(gb_text: str) -> dict | None:
    """Extract CDS start/end positions and translation from GenBank features.

    Returns dict with 'start' (1-based mRNA position), 'end' (1-based inclusive),
    and 'protein_start' (first 30 aa of translation), or None if no CDS found.
    """
    record = SeqIO.read(StringIO(gb_text), "genbank")
    for feature in record.features:
        if feature.type == "CDS":
            cds_start = int(feature.location.start) + 1
            cds_end = int(feature.location.end)
            translation = feature.qualifiers.get("translation", [""])[0]
            return {
                "start": cds_start,
                "end": cds_end,
                "protein_start": translation[:30],
            }
    return None


def parse_exon_features(gb_text: str) -> list[tuple[int, int]]:
    """Extract sorted exon (start, end) positions from GenBank features.

    Returns 1-based inclusive positions on the mRNA.
    """
    record = SeqIO.read(StringIO(gb_text), "genbank")
    exons = []
    for feature in record.features:
        if feature.type == "exon":
            start = int(feature.location.start) + 1  # convert 0-based to 1-based
            end = int(feature.location.end)           # 1-based inclusive
            exons.append((start, end))
    exons.sort(key=lambda x: x[0])
    return exons


def compute_junctions(exons: list[tuple[int, int]]) -> list[dict]:
    """Compute exon-exon junction positions from consecutive exons.

    The junction mrna_position is the position of the last nucleotide
    of the upstream exon (equivalently, one before the first nt of the next exon).
    """
    junctions = []
    for i in range(len(exons) - 1):
        exon_before_end = exons[i][1]
        junctions.append({
            "junction_id": i + 1,
            "mrna_position": exon_before_end,
            "exon_before": i + 1,
            "exon_after": i + 2,
        })
    return junctions


def main():
    api_key = snakemake.params.api_key
    transcript_id_file = snakemake.input.transcript_id
    output_transcript = snakemake.output.transcript
    output_junctions = snakemake.output.junctions

    transcript_id = Path(transcript_id_file).read_text().strip()
    log.info(f"Fetching transcript data for {transcript_id}")

    gb_text = fetch_genbank(transcript_id, api_key)

    record = SeqIO.read(StringIO(gb_text), "genbank")
    sequence = str(record.seq)
    accession = record.id
    organism = record.annotations.get("organism", "")

    log.info(
        f"Transcript: {accession}, length: {len(sequence)}, organism: {organism}"
    )

    # Parse CDS feature
    cds_info = parse_cds_feature(gb_text)
    if cds_info:
        log.info(
            f"CDS: mRNA positions {cds_info['start']}-{cds_info['end']}, "
            f"protein starts with {cds_info['protein_start']}"
        )
    else:
        log.warning("No CDS feature found in GenBank record")

    # Write transcript TSV
    with open(output_transcript, "w") as f:
        f.write("field\tvalue\n")
        f.write(f"accession\t{accession}\n")
        f.write(f"sequence_length\t{len(sequence)}\n")
        f.write(f"organism\t{organism}\n")
        if cds_info:
            f.write(f"cds_start\t{cds_info['start']}\n")
            f.write(f"cds_end\t{cds_info['end']}\n")
            f.write(f"protein_start\t{cds_info['protein_start']}\n")
        f.write(f"sequence\t{sequence}\n")
    log.info(f"Wrote transcript to {output_transcript}")

    # Parse exon features and compute junctions
    exons = parse_exon_features(gb_text)
    log.info(f"Found {len(exons)} exons")

    junctions = compute_junctions(exons)
    log.info(f"Computed {len(junctions)} exon-exon junctions")

    # Write junctions TSV
    with open(output_junctions, "w") as f:
        f.write("junction_id\tmrna_position\texon_before\texon_after\n")
        for j in junctions:
            f.write(
                f"{j['junction_id']}\t{j['mrna_position']}\t"
                f"{j['exon_before']}\t{j['exon_after']}\n"
            )
    log.info(f"Wrote junctions to {output_junctions}")


if __name__ == "__main__":
    main()
