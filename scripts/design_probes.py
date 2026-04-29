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


def design_probes(
    mutation_positions: list[int],
    mutation_labels: list[str],
    sequence: str,
    junction_positions: list[int],
    target_len: int,
    min_len: int,
    max_len: int,
    junction_buffer: int,
    scan_offset: int,
    gene: str,
) -> list[dict]:
    forbidden = build_forbidden_zones(junction_positions, junction_buffer)
    log.info(f"Forbidden zones: {forbidden}")

    seq_len = len(sequence)
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

                covered = [
                    idx
                    for idx in uncovered
                    if probe_start <= mutation_positions[idx] <= probe_end
                ]

                overlaps_existing = any(
                    intervals_overlap(probe_start, probe_end, p["start"], p["end"])
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
                f"at position {target_pos} — skipping"
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

        probes.append({
            "probe_name": f"{gene}-{probe_counter}",
            "start": best_probe["start"],
            "end": best_probe["end"],
            "length": best_probe["length"],
            "target_sequence": target_seq,
            "fwd_oligo": fwd_oligo,
            "rev_oligo": rev_oligo,
            "mutations_covered": ";".join(covered_mutations),
            "mutation_positions": ";".join(str(p) for p in covered_positions),
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

    mutations_df = pd.read_csv(snakemake.input.mutations, sep="\t")
    transcript_df = pd.read_csv(snakemake.input.transcript, sep="\t")
    junctions_df = pd.read_csv(snakemake.input.junctions, sep="\t")

    sequence = transcript_df.loc[transcript_df["field"] == "sequence", "value"].values[0]
    junction_positions = junctions_df["mrna_position"].tolist()

    mutation_positions = mutations_df["mrna_position"].tolist()
    mutation_labels = mutations_df["protein_change"].tolist()

    log.info(
        f"Designing probes for {gene}: {len(mutation_positions)} mutations, "
        f"{len(junction_positions)} junctions, sequence length {len(sequence)}"
    )

    probes = design_probes(
        mutation_positions=mutation_positions,
        mutation_labels=mutation_labels,
        sequence=sequence,
        junction_positions=junction_positions,
        target_len=target_len,
        min_len=min_len,
        max_len=max_len,
        junction_buffer=junction_buffer,
        scan_offset=scan_offset,
        gene=gene,
    )

    log.info(f"Designed {len(probes)} probes")

    covered = set()
    for p in probes:
        for pos_str in p["mutation_positions"].split(";"):
            covered.add(int(pos_str))
    total = set(mutation_positions)
    log.info(
        f"Coverage: {len(covered)}/{len(total)} mutations "
        f"({100 * len(covered) / len(total):.1f}%)"
    )

    probes_df = pd.DataFrame(probes)
    probes_df.to_csv(snakemake.output[0], sep="\t", index=False)
    log.info(f"Wrote probes to {snakemake.output[0]}")


if __name__ == "__main__":
    main()
