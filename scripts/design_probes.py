"""Design antisense probes covering missense mutations on a transcript.

Greedy algorithm: sweep left to right, placing probes that maximize
mutation coverage while avoiding exon-exon junction buffer zones.
"""

import logging
import sys

import pandas as pd

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent intervals."""
    if not intervals:
        return []
    sorted_iv = sorted(intervals)
    merged = [sorted_iv[0]]
    for start, end in sorted_iv[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def overlaps_any(start: int, end: int, forbidden: list[tuple[int, int]]) -> bool:
    return any(intervals_overlap(start, end, fs, fe) for fs, fe in forbidden)


def build_forbidden_zones(junctions: list[int], buffer: int) -> list[tuple[int, int]]:
    """Build forbidden zones around each exon-exon junction.

    A junction at position P means the last nt of exon i is at P and
    the first nt of exon i+1 is at P+1. We exclude [P-buffer+1, P+buffer].
    """
    zones = []
    for pos in junctions:
        zones.append((pos - buffer + 1, pos + buffer))
    return merge_intervals(zones)


def compute_probe_exon(mrna_start: int, mrna_end: int,
                       junction_positions: list[int],
                       cds_mrna_start: int, cds_mrna_end: int) -> str:
    """Determine the exon label for a probe based on its midpoint."""
    mid = (mrna_start + mrna_end) // 2
    if mid < cds_mrna_start:
        return "5UTR"
    if mid > cds_mrna_end:
        return "3UTR"
    for i, jpos in enumerate(junction_positions):
        if mid <= jpos:
            return str(i + 1)
    return str(len(junction_positions) + 1)


def parse_region_config(region_str: str) -> list[str]:
    """Parse region config string into list of region tokens.

    Examples:
        "CDS"       -> ["CDS"]
        "All"       -> ["All"]
        "5"         -> ["5"]
        "3-8"       -> ["3", "4", "5", "6", "7", "8"]
        "1-4,8-10"  -> ["1", "2", "3", "4", "8", "9", "10"]
        "5UTR,3-5"  -> ["5UTR", "3", "4", "5"]
    """
    tokens = []
    for token in region_str.split(","):
        token = token.strip()
        if not token:
            continue
        if token in ("5UTR", "3UTR", "CDS", "All"):
            tokens.append(token)
        elif "-" in token:
            parts = token.split("-")
            start, end = int(parts[0]), int(parts[1])
            for n in range(start, end + 1):
                tokens.append(str(n))
        else:
            tokens.append(token)
    return tokens


def build_exon_boundaries(junction_positions: list[int], seq_len: int) -> list[tuple[int, int]]:
    """Build exon boundary list [(start, end), ...] from junction positions."""
    bounds = []
    prev_end = 0
    for jpos in junction_positions:
        bounds.append((prev_end + 1, jpos))
        prev_end = jpos
    bounds.append((prev_end + 1, seq_len))
    return bounds


def mrna_ranges_for_regions(region_tokens: list[str],
                            junction_positions: list[int],
                            seq_len: int,
                            cds_mrna_start: int,
                            cds_mrna_end: int) -> list[tuple[int, int]]:
    """Convert region tokens to mRNA position ranges.

    Returns list of (start, end) 1-based inclusive ranges on the mRNA.
    Raises ValueError for out-of-range exon numbers.
    """
    n_exons = len(junction_positions) + 1
    exon_bounds = build_exon_boundaries(junction_positions, seq_len)

    ranges = []
    for token in set(region_tokens):
        if token == "5UTR":
            if cds_mrna_start > 1:
                ranges.append((1, cds_mrna_start - 1))
        elif token == "3UTR":
            if cds_mrna_end < seq_len:
                ranges.append((cds_mrna_end + 1, seq_len))
        elif token == "CDS":
            ranges.append((cds_mrna_start, cds_mrna_end))
        elif token == "All":
            ranges.append((1, seq_len))
        else:
            exon_num = int(token)
            if exon_num < 1 or exon_num > n_exons:
                raise ValueError(
                    f"Exon {exon_num} is out of range. "
                    f"Transcript has {n_exons} exons (1-{n_exons}). "
                    f"Please update the 'region' setting in config.yaml."
                )
            ranges.append(exon_bounds[exon_num - 1])

    return merge_intervals(ranges)


def design_probes(
    mutation_positions: list[int],
    mutation_labels: list[str],
    mutation_cdna_changes: list[str],
    sequence: str,
    junction_positions: list[int],
    target_len: int,
    min_len: int,
    max_len: int,
    junction_buffer: int,
    scan_offset: int,
    gene: str,
    cds_mrna_start: int,
    cds_mrna_end: int,
    region_tokens: list[str] | None = None,
) -> list[dict]:
    forbidden = build_forbidden_zones(junction_positions, junction_buffer)
    log.info(f"Forbidden zones: {forbidden}")

    seq_len = len(sequence)

    # Build allowed zones from region config
    if region_tokens is not None:
        allowed = mrna_ranges_for_regions(
            region_tokens, junction_positions, seq_len,
            cds_mrna_start, cds_mrna_end,
        )
        log.info(f"Region allowed zones: {allowed}")
    else:
        allowed = [(1, seq_len)]
    uncovered = set(range(len(mutation_positions)))
    probes = []
    probe_counter = 0

    while uncovered:
        first_idx = min(uncovered, key=lambda i: mutation_positions[i])
        target_pos = mutation_positions[first_idx]

        best_probe = None
        best_score = -1

        # Try lengths: target first, then expanding outward
        lengths = [target_len]
        for delta in range(1, max_len - target_len + 1):
            if target_len - delta >= min_len:
                lengths.append(target_len - delta)
            if target_len + delta <= max_len:
                lengths.append(target_len + delta)

        for length in lengths:
            for probe_start in range(
                max(scan_offset + 1, target_pos - length + 1),
                target_pos + 1,
            ):
                probe_end = probe_start + length - 1  # 1-based inclusive

                if probe_end > seq_len:
                    continue

                if overlaps_any(probe_start, probe_end, forbidden):
                    continue

                # Check if probe falls within allowed regions
                if not any(
                    probe_start >= a_start and probe_end <= a_end
                    for a_start, a_end in allowed
                ):
                    continue

                covered = [
                    idx
                    for idx in uncovered
                    if probe_start <= mutation_positions[idx] <= probe_end
                ]

                overlaps_existing = any(
                    intervals_overlap(probe_start, probe_end, p["mrna_start"], p["mrna_end"])
                    for p in probes
                )

                if covered:
                    min_pos = min(mutation_positions[i] for i in covered)
                    max_pos = max(mutation_positions[i] for i in covered)
                    mutation_mid = (min_pos + max_pos) / 2
                else:
                    mutation_mid = target_pos

                center = (probe_start + probe_end) / 2
                centering = 1.0 - abs(mutation_mid - center) / (length / 2)

                score = len(covered) * 100 + centering * 10
                if overlaps_existing:
                    score -= 50

                if score > best_score:
                    best_score = score
                    best_probe = {
                        "start": probe_start,
                        "end": probe_end,
                        "length": length,
                        "covered_indices": covered,
                    }

        if best_probe is None:
            log.warning(
                f"Cannot place probe for mutation {mutation_labels[first_idx]} "
                f"at mRNA position {target_pos} — skipping"
            )
            uncovered.discard(first_idx)
            continue

        probe_counter += 1
        covered_idxs = best_probe["covered_indices"]

        s = best_probe["start"] - 1  # 0-based
        e = best_probe["end"]         # exclusive upper bound
        target_seq = sequence[s:e]

        fwd_oligo = "AAAC" + reverse_complement(target_seq)
        rev_oligo = "AAAA" + target_seq

        covered_mutations = [mutation_labels[i] for i in covered_idxs]
        covered_positions = [mutation_positions[i] for i in covered_idxs]
        covered_cdna = [mutation_cdna_changes[i] for i in covered_idxs]

        mrna_start = best_probe["start"]
        mrna_end = best_probe["end"]
        cds_start = mrna_start - cds_mrna_start + 1
        cds_end = mrna_end - cds_mrna_start + 1

        exon_label = compute_probe_exon(
            mrna_start, mrna_end, junction_positions,
            cds_mrna_start, cds_mrna_end,
        )

        probes.append({
            "probe_name": f"{gene}-{exon_label}-{probe_counter}",
            "mrna_start": mrna_start,
            "mrna_end": mrna_end,
            "cds_start": cds_start,
            "cds_end": cds_end,
            "length": best_probe["length"],
            "target_sequence": target_seq,
            "fwd_oligo": fwd_oligo,
            "rev_oligo": rev_oligo,
            "mutations_covered": ",".join(covered_mutations),
            "cdna_changes_covered": ",".join(covered_cdna),
            "mutation_positions": ",".join(str(p) for p in covered_positions),
            "exon": exon_label,
        })

        uncovered -= set(covered_idxs)

    return probes


def main():
    gene = snakemake.params.gene
    target_len = snakemake.params.target_length
    min_len = snakemake.params.min_length
    max_len = snakemake.params.max_length
    junction_buffer = snakemake.params.junction_buffer
    scan_offset = snakemake.params.scan_offset
    region_config = snakemake.params.region

    mutations_df = pd.read_csv(snakemake.input.mutations, sep="\t")
    transcript_df = pd.read_csv(snakemake.input.transcript, sep="\t")
    junctions_df = pd.read_csv(snakemake.input.junctions, sep="\t")

    sequence = transcript_df.loc[transcript_df["field"] == "sequence", "value"].values[0]
    junction_positions = junctions_df["mrna_position"].tolist()

    # Read CDS start from transcript.tsv
    cds_start_row = transcript_df.loc[transcript_df["field"] == "cds_start", "value"]
    if cds_start_row.empty:
        log.error("No cds_start found in transcript.tsv — re-run fetch_transcript")
        sys.exit(1)
    cds_mrna_start = int(cds_start_row.values[0])
    log.info(f"CDS starts at mRNA position {cds_mrna_start}")

    # Read CDS end from transcript.tsv
    cds_end_row = transcript_df.loc[transcript_df["field"] == "cds_end", "value"]
    if cds_end_row.empty:
        log.error("No cds_end found in transcript.tsv — re-run fetch_transcript")
        sys.exit(1)
    cds_mrna_end = int(cds_end_row.values[0])
    log.info(f"CDS ends at mRNA position {cds_mrna_end}")

    # Convert cdna_position to mrna_position, handling UTR
    mutations_with_pos = mutations_df[mutations_df["cdna_position"].notna()].copy()
    mutations_with_pos["cdna_position"] = mutations_with_pos["cdna_position"].astype(int)

    def cdna_to_mrna(row):
        cdna_type = row.get("cdna_type")
        cdna_pos = int(row["cdna_position"])
        if cdna_type == "5UTR":
            return cds_mrna_start - cdna_pos
        elif cdna_type == "3UTR":
            return cds_mrna_end + cdna_pos
        else:
            return cdna_pos + cds_mrna_start - 1

    mutations_with_pos["mrna_position"] = mutations_with_pos.apply(cdna_to_mrna, axis=1)

    n_total = len(mutations_df)
    n_without_cdna = n_total - len(mutations_with_pos)
    if n_without_cdna > 0:
        log.info(
            f"{n_without_cdna}/{n_total} mutations have no valid cDNA position — "
            f"excluded from probe design"
        )

    # Filter mutations by region config
    if region_config and region_config != "All":
        region_tokens = parse_region_config(region_config)
        allowed_ranges = mrna_ranges_for_regions(
            region_tokens, junction_positions, len(sequence),
            cds_mrna_start, cds_mrna_end,
        )
        n_before = len(mutations_with_pos)

        def mutation_in_region(mrna_pos):
            return any(s <= mrna_pos <= e for s, e in allowed_ranges)

        mutations_with_pos = mutations_with_pos[
            mutations_with_pos["mrna_position"].apply(mutation_in_region)
        ]
        n_after = len(mutations_with_pos)
        if n_after < n_before:
            log.info(
                f"Region filter '{region_config}': kept {n_after}/{n_before} mutations"
            )
    else:
        region_tokens = None

    mutations_with_pos["protein_change"] = mutations_with_pos["protein_change"].fillna("")
    mutations_with_pos["cdna_change"] = mutations_with_pos["cdna_change"].fillna("")

    mutation_positions = mutations_with_pos["mrna_position"].tolist()
    mutation_labels = mutations_with_pos["protein_change"].tolist()
    mutation_cdna_changes = mutations_with_pos["cdna_change"].tolist()

    log.info(
        f"Designing probes for {gene}: {len(mutation_positions)} mutations, "
        f"{len(junction_positions)} junctions, sequence length {len(sequence)}"
    )

    probes = design_probes(
        mutation_positions=mutation_positions,
        mutation_labels=mutation_labels,
        mutation_cdna_changes=mutation_cdna_changes,
        sequence=sequence,
        junction_positions=junction_positions,
        target_len=target_len,
        min_len=min_len,
        max_len=max_len,
        junction_buffer=junction_buffer,
        scan_offset=scan_offset,
        gene=gene,
        cds_mrna_start=cds_mrna_start,
        cds_mrna_end=cds_mrna_end,
        region_tokens=region_tokens,
    )

    log.info(f"Designed {len(probes)} probes")

    covered = set()
    for p in probes:
        for pos_str in p["mutation_positions"].split(","):
            covered.add(int(pos_str))
    total = set(mutation_positions)
    if total:
        log.info(
            f"Coverage: {len(covered)}/{len(total)} mutations "
            f"({100 * len(covered) / len(total):.1f}%)"
        )
    else:
        log.info("No mutations with valid positions — no probes designed")

    probes_df = pd.DataFrame(probes)
    probes_df.to_csv(snakemake.output[0], sep="\t", index=False)
    log.info(f"Wrote probes to {snakemake.output[0]}")


if __name__ == "__main__":
    main()
