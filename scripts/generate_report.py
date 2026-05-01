"""Generate interactive HTML report with Plotly visualization."""

import json
import logging
from datetime import date

import pandas as pd
import plotly.graph_objects as go
from jinja2 import Template

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)

SIG_COLORS = {
    "Pathogenic": "#e74c3c",
    "Pathogenic/Likely pathogenic": "#e74c3c",
    "Likely pathogenic": "#e67e22",
    "Uncertain significance": "#f39c12",
    "Likely benign": "#27ae60",
    "Benign": "#2ecc71",
    "Benign/Likely benign": "#2ecc71",
}

SIG_CLASS = {
    "Pathogenic": "p",
    "Pathogenic/Likely pathogenic": "p",
    "Likely pathogenic": "lp",
    "Uncertain significance": "v",
    "Likely benign": "lb",
    "Benign": "b",
    "Benign/Likely benign": "b",
}


def load_data(mutations_path, transcript_path, junctions_path, probes_path):
    mutations_df = pd.read_csv(mutations_path, sep="\t")
    transcript_df = pd.read_csv(transcript_path, sep="\t")
    junctions_df = pd.read_csv(junctions_path, sep="\t")
    probes_df = pd.read_csv(probes_path, sep="\t")

    def get_field(name):
        row = transcript_df.loc[transcript_df["field"] == name, "value"]
        return row.values[0] if not row.empty else None

    sequence = get_field("sequence")
    accession = get_field("accession")
    cds_start = int(get_field("cds_start")) if get_field("cds_start") else None
    cds_end = int(get_field("cds_end")) if get_field("cds_end") else None

    return mutations_df, sequence, accession, cds_start, cds_end, junctions_df, probes_df


def compute_exon_number(mrna_pos, junctions_df, cds_start=None, cds_end=None):
    """Determine which exon a mRNA position falls in.

    Returns numeric exon string, or '5UTR'/'3UTR' for UTR positions.
    """
    if pd.isna(mrna_pos):
        return ""
    pos = int(mrna_pos)
    if cds_start and pos < cds_start:
        return "5UTR"
    if cds_end and pos > cds_end:
        return "3UTR"
    junction_positions = junctions_df["mrna_position"].tolist()
    for i, jpos in enumerate(junction_positions):
        if pos <= jpos:
            return str(i + 1)
    return str(len(junction_positions) + 1)


def _parse_region_config(region_str):
    """Parse region config string into list of region tokens."""
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


def _mrna_ranges_for_regions(region_tokens, junction_positions, seq_len,
                              cds_start, cds_end):
    """Convert region tokens to mRNA position ranges."""
    n_exons = len(junction_positions) + 1
    exon_bounds = []
    prev_end = 0
    for jpos in junction_positions:
        exon_bounds.append((prev_end + 1, jpos))
        prev_end = jpos
    exon_bounds.append((prev_end + 1, seq_len))

    ranges = []
    for token in set(region_tokens):
        if token == "5UTR":
            if cds_start > 1:
                ranges.append((1, cds_start - 1))
        elif token == "3UTR":
            if cds_end < seq_len:
                ranges.append((cds_end + 1, seq_len))
        elif token == "CDS":
            ranges.append((cds_start, cds_end))
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

    # Merge overlapping ranges
    if not ranges:
        return []
    sorted_iv = sorted(ranges)
    merged = [sorted_iv[0]]
    for start, end in sorted_iv[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def build_enriched_mutations(mutations_df, cds_start, cds_end, junctions_df,
                             probes_df):
    """Add computed columns to mutations: mrna_position, exon, probe_name, id."""
    df = mutations_df.copy()

    # Compute mRNA position from cdna_position, handling UTR
    if cds_start and cds_end and "cdna_position" in df.columns:
        has_pos = df["cdna_position"].notna()
        pos_df = df[has_pos].copy()
        pos_df["cdna_position"] = pos_df["cdna_position"].astype(int)

        def _cdna_to_mrna(row):
            cdna_type = row.get("cdna_type")
            cdna_pos = int(row["cdna_position"])
            if cdna_type == "5UTR":
                return cds_start - cdna_pos
            elif cdna_type == "3UTR":
                return cds_end + cdna_pos
            else:
                return cdna_pos + cds_start - 1

        df.loc[has_pos, "mrna_position"] = pos_df.apply(_cdna_to_mrna, axis=1)
        df.loc[~has_pos, "mrna_position"] = None
    elif cds_start and "cdna_position" in df.columns:
        has_pos = df["cdna_position"].notna()
        df.loc[has_pos, "mrna_position"] = (
            df.loc[has_pos, "cdna_position"].astype(int) + cds_start - 1
        ).astype(int)
        df.loc[~has_pos, "mrna_position"] = None
    elif "mrna_position" not in df.columns:
        df["mrna_position"] = None

    # HGVS ID: NM_007055.4(POLR3A):c.251G>A
    df["id"] = df["transcript_id"] + "(" + df["gene"] + "):" + df["cdna_change"]

    # Exon number
    df["exon"] = df["mrna_position"].apply(
        lambda p: compute_exon_number(p, junctions_df, cds_start, cds_end)
    )

    # Probe name: which probe covers this mutation
    probe_lookup = {}
    for _, prow in probes_df.iterrows():
        for pos_str in str(prow["mutation_positions"]).split(","):
            if pos_str.strip():
                pos = int(pos_str.strip())
                probe_lookup.setdefault(pos, []).append(prow["probe_name"])

    df["probe"] = df["mrna_position"].apply(
        lambda p: ", ".join(sorted(set(probe_lookup.get(int(p), [])))) if pd.notna(p) and int(p) in probe_lookup else ""
    )

    return df


def build_plotly_figure(mutations_df, accession, junctions_df, probes_df, seq_len):
    fig = go.Figure()

    # Exon blocks
    exons = []
    prev_end = 0
    for _, jrow in junctions_df.iterrows():
        exons.append((prev_end + 1, jrow["mrna_position"]))
        prev_end = jrow["mrna_position"]
    exons.append((prev_end + 1, seq_len))

    n_exons = len(exons)
    for i, (start, end) in enumerate(exons):
        color = f"hsl({int(210 + 40 * (i / n_exons))}, 60%, 75%)"
        fig.add_trace(
            go.Bar(
                x=[end - start + 1], y=[1], base=start - 1,
                orientation="h", marker_color=color,
                hovertemplate=f"Exon {i + 1}<br>Position: {start}-{end}<extra></extra>",
                showlegend=False,
            )
        )

    # Mutations (only those with valid mrna_position)
    mut_valid = mutations_df[mutations_df["mrna_position"].notna()]
    mut_colors = [SIG_COLORS.get(s, "#95a5a6") for s in mut_valid["clinical_significance"]]
    hover_text = [
        f"<b>{row['protein_change']}</b><br>"
        f"cDNA: {row['cdna_change']}<br>"
        f"Significance: {row['clinical_significance']}<br>"
        f"Accession: {row['clinvar_accession']}<br>"
        f"mRNA Position: {int(row['mrna_position'])}<br>"
        f"CDS Position: {int(row['cdna_position']) if pd.notna(row.get('cdna_position')) else ''}"
        for _, row in mut_valid.iterrows()
    ]
    fig.add_trace(
        go.Scatter(
            x=mut_valid["mrna_position"].astype(int).tolist(),
            y=[1.8] * len(mut_valid),
            mode="markers",
            marker=dict(size=6, color=mut_colors, line=dict(width=0.5, color="white")),
            text=hover_text, hoverinfo="text", showlegend=False,
        )
    )

    # Probes
    for _, prow in probes_df.iterrows():
        muts = prow["mutations_covered"]
        fig.add_trace(
            go.Bar(
                x=[prow["mrna_end"] - prow["mrna_start"] + 1],
                y=[2.6], base=prow["mrna_start"] - 1,
                orientation="h",
                marker_color="rgba(52, 152, 219, 0.5)",
                marker_line=dict(color="rgba(41, 128, 185, 0.8)", width=1),
                hovertemplate=(
                    f"<b>{prow['probe_name']}</b><br>"
                    f"mRNA: {prow['mrna_start']}-{prow['mrna_end']} ({prow['length']} nt)<br>"
                    f"CDS: {prow['cds_start']}-{prow['cds_end']}<br>"
                    f"Mutations: {muts}<extra></extra>"
                ),
                showlegend=False,
            )
        )

    fig.update_layout(
        barmode="overlay",
        xaxis=dict(title="mRNA Position (nt)", rangeslider=dict(visible=True)),
        yaxis=dict(tickvals=[1, 1.8, 2.6], ticktext=["Exons", "Mutations", "Probes"], range=[0.3, 3.2]),
        height=400, margin=dict(l=80, r=30, t=30, b=60), plot_bgcolor="white",
    )
    return fig


def build_annotated_sequence(sequence, cds_start, cds_end, junctions_df, mutations_df, probes_df):
    """Generate HTML for the annotated mRNA sequence using CSS classes.
    Returns (visual_html, copy_html) — two versions of the annotated sequence.
    """
    seq_len = len(sequence)
    nts_per_line = 100
    copy_nts_per_line = 60

    # Build exon boundaries (1-based inclusive)
    exon_bounds = []
    prev_end = 0
    for _, jrow in junctions_df.iterrows():
        exon_bounds.append((prev_end + 1, int(jrow["mrna_position"])))
        prev_end = int(jrow["mrna_position"])
    exon_bounds.append((prev_end + 1, seq_len))

    # Map mutation mRNA position -> list of (sig_class, color_hex, tooltip)
    mut_map = {}
    for _, mrow in mutations_df.iterrows():
        if pd.notna(mrow.get("mrna_position")):
            pos = int(mrow["mrna_position"])
            sig = mrow["clinical_significance"]
            sig_cls = SIG_CLASS.get(sig, "o")
            color = SIG_COLORS.get(sig, "#e74c3c")
            tip = f"{mrow['id']} | {sig} | mRNA pos: {pos}"
            mut_map.setdefault(pos, []).append((sig_cls, color, tip))

    # Probe colors: alternating
    PROBE_COLORS = ["#dee6ef", "#ffffd7"]
    probe_color_map = {}
    for i, (_, prow) in enumerate(probes_df.iterrows()):
        probe_color_map[prow["probe_name"]] = PROBE_COLORS[i % 2]

    def blend_colors(colors):
        if len(colors) == 1:
            return colors[0]
        r = sum(int(c[1:3], 16) for c in colors) // len(colors)
        g = sum(int(c[3:5], 16) for c in colors) // len(colors)
        b = sum(int(c[5:7], 16) for c in colors) // len(colors)
        return f"#{r:02x}{g:02x}{b:02x}"

    # Map mRNA position -> background color from probe(s)
    pos_probes = {}
    for _, prow in probes_df.iterrows():
        name = prow["probe_name"]
        for p in range(int(prow["mrna_start"]), int(prow["mrna_end"]) + 1):
            pos_probes.setdefault(p, []).append(name)

    probe_bg = {}
    for p, names in pos_probes.items():
        unique_colors = list(set(probe_color_map[n] for n in names))
        probe_bg[p] = blend_colors(unique_colors)

    # Probe info for label placement
    probe_info = []
    for _, prow in probes_df.iterrows():
        probe_info.append({
            "name": prow["probe_name"],
            "start": int(prow["mrna_start"]),
            "end": int(prow["mrna_end"]),
        })

    def get_labels_for_line(line_start, line_end):
        """Return [(midpoint_position, name), ...] for probes showing labels."""
        labels = []
        for pi in probe_info:
            if pi["end"] < line_start or pi["start"] > line_end:
                continue
            overlap_start = max(pi["start"], line_start)
            overlap_end = min(pi["end"], line_end)
            nts_on_line = overlap_end - overlap_start + 1
            total_probe_nts = pi["end"] - pi["start"] + 1
            mid = (pi["start"] + pi["end"]) // 2

            if pi["start"] >= line_start and pi["end"] <= line_end:
                label_pos = mid
            else:
                nts_other = total_probe_nts - nts_on_line
                if nts_on_line > nts_other:
                    label_pos = (overlap_start + overlap_end) // 2
                elif nts_on_line < nts_other:
                    continue
                else:
                    if line_start <= pi["start"] or line_start < mid:
                        label_pos = (overlap_start + overlap_end) // 2
                    else:
                        continue
            labels.append((label_pos, pi["name"]))
        return labels

    def render_label_line(line_start, line_end, labels):
        """Render a label line using invisible dashes for character-level alignment."""
        n_nts = line_end - line_start + 1
        n_groups = (n_nts + 9) // 10

        # Build character array: all positions are invisible dashes (including separators)
        chars = []
        for g in range(n_groups):
            group_size = min(10, n_nts - g * 10)
            for _ in range(group_size):
                chars.append(("p", "-"))
            if g < n_groups - 1:
                chars.append(("p", "-"))

        total_visual = len(chars)

        # Overlay labels at their center positions
        for midpoint, name in labels:
            char_offset = midpoint - line_start
            visual_center = char_offset + char_offset // 10
            label_start = visual_center - len(name) // 2

            for i, ch in enumerate(name):
                vpos = label_start + i
                if 0 <= vpos < total_visual and chars[vpos][0] == "p":
                    chars[vpos] = ("l", ch)

        # Render by grouping consecutive same-type characters
        parts = []
        i = 0
        while i < len(chars):
            typ = chars[i][0]
            j = i
            while j < len(chars) and chars[j][0] == typ:
                j += 1
            text = "".join(chars[k][1] for k in range(i, j))
            cls = "ll" if typ == "l" else "lp"
            parts.append(f'<span class="{cls}">{text}</span>')
            i = j

        return "".join(parts)

    def render_nucleotide_line(pos, line_end, for_copy=False):
        """Render one line of nucleotides. Visual uses CSS classes, copy uses inline styles."""
        parts = []
        for group_start in range(pos, line_end + 1, 10):
            group_end = min(group_start + 9, line_end)
            for p in range(group_start, group_end + 1):
                nt = sequence[p - 1]
                is_utr = (cds_start and p < cds_start) or (cds_end and p > cds_end)

                if for_copy:
                    # Inline styles for rich text paste compatibility
                    styles = "font-family:'Courier New',monospace;font-size:12px;"
                    if p in probe_bg:
                        styles += f"background:{'#e0e0e0' if is_utr else probe_bg[p]};"
                    elif is_utr:
                        styles += "background:#e0e0e0;"
                    if p in mut_map:
                        styles += f"color:{mut_map[p][0][1]};font-weight:bold;"
                    parts.append(f'<span style="{styles}">{nt}</span>')
                else:
                    # CSS classes for visual HTML
                    classes = ["s"]
                    bg_style = ""
                    if p in probe_bg:
                        if is_utr:
                            classes.append("su")
                        else:
                            bg_style = f' style="background:{probe_bg[p]}"'
                    elif is_utr:
                        classes.append("su")
                    if p in mut_map:
                        classes.append("sm")
                        classes.append(f'sm-{mut_map[p][0][0]}')

                    cls_str = " ".join(classes)
                    if p in mut_map:
                        tips = "; ".join(t for _, _, t in mut_map[p])
                        tips_escaped = tips.replace("&", "&amp;").replace('"', "&quot;").replace("<", "&lt;").replace(">", "&gt;")
                        parts.append(f'<span class="mutation-tip {cls_str}"{bg_style} data-tip="{tips_escaped}">{nt}</span>')
                    else:
                        parts.append(f'<span class="{cls_str}"{bg_style}>{nt}</span>')

            if group_end < line_end:
                if for_copy:
                    parts.append('<span style="font-family:\'Courier New\',monospace;"> </span>')
                else:
                    parts.append('<span class="sg"> </span>')
        return parts

    # Build visual HTML
    visual_parts = []
    for exon_idx, (ex_start, ex_end) in enumerate(exon_bounds, 1):
        visual_parts.append(
            f'<div class="se">'
            f'<div class="sh">'
            f'Exon {exon_idx} (positions {ex_start}-{ex_end})</div>'
        )
        pos = ex_start
        while pos <= ex_end:
            line_end = min(pos + nts_per_line - 1, ex_end)

            # Probe labels
            labels = get_labels_for_line(pos, line_end)
            if labels:
                label_html = render_label_line(pos, line_end, labels)
                visual_parts.append(f'<div class="sl">{label_html}</div>')
            else:
                visual_parts.append('<div class="sp"></div>')

            # Nucleotides
            nt_parts = render_nucleotide_line(pos, line_end, for_copy=False)
            visual_parts.append('<div class="sn">' + "".join(nt_parts) + "</div>")
            pos = line_end + 1

        visual_parts.append("</div>")

    # Build copy HTML
    copy_parts = []
    for exon_idx, (ex_start, ex_end) in enumerate(exon_bounds, 1):
        exon_probes = []
        for pi in probe_info:
            if pi["end"] >= ex_start and pi["start"] <= ex_end:
                exon_probes.append(pi["name"])
        header = f"Exon {exon_idx} (positions {ex_start}-{ex_end}"
        if exon_probes:
            probe_nums = []
            for pn in exon_probes:
                parts = pn.split("-")
                # New format: {gene}-{exon}-{N} -> last part is the number
                last_part = parts[-1]
                try:
                    probe_nums.append(int(last_part))
                except ValueError:
                    probe_nums.append(pn)
            if all(isinstance(n, int) for n in probe_nums):
                probe_nums.sort()
                if probe_nums[-1] - probe_nums[0] == len(probe_nums) - 1:
                    header += f", probes {probe_nums[0]} to {probe_nums[-1]}"
                else:
                    header += f", probes {', '.join(str(n) for n in probe_nums)}"
        header += ")"
        copy_parts.append(
            f'<div class="se">'
            f'<div class="sh">{header}</div>'
        )
        pos = ex_start
        while pos <= ex_end:
            line_end = min(pos + copy_nts_per_line - 1, ex_end)
            nt_parts = render_nucleotide_line(pos, line_end, for_copy=True)
            copy_parts.append('<div class="sn">' + "".join(nt_parts) + "</div>")
            pos = line_end + 1
        copy_parts.append("</div>")
        if exon_idx < len(exon_bounds):
            copy_parts.append('<p style="margin:6px 0;">&nbsp;</p>')

    return "\n".join(visual_parts), "\n".join(copy_parts)


def write_mutations_file(mutations_df, output_path):
    """Write the detailed mutations TSV file."""
    cols = ["id", "mrna_position", "cdna_position", "clinical_significance",
            "cdna_change", "protein_change", "clinvar_accession",
            "chrom", "pos_start", "exon", "probe"]
    out_df = mutations_df[[c for c in cols if c in mutations_df.columns]].copy()
    out_df = out_df.rename(columns={
        "mrna_position": "mRNA Position",
        "cdna_position": "CDS Position",
        "clinical_significance": "Significance",
        "cdna_change": "cDNA Change",
        "protein_change": "Protein Change",
        "clinvar_accession": "ClinVar Accession",
        "chrom": "Chrom",
        "pos_start": "Genomic Pos",
        "exon": "Exon",
        "probe": "Probe",
        "id": "ID",
    })
    out_df.to_csv(output_path, sep="\t", index=False)


HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{{ gene }} Probe Design Report</title>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<style>
body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px 40px; color: #2c3e50; }
h1 { border-bottom: 2px solid #3498db; padding-bottom: 8px; }
h2 { color: #2c3e50; margin-top: 30px; }
.section { margin-bottom: 40px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { border: 1px solid #ddd; padding: 6px 10px; text-align: left; white-space: nowrap; }
th { background-color: #2c3e50; color: white; cursor: pointer; user-select: none; }
th:hover { background-color: #34495e; }
tr:nth-child(even) { background-color: #f9f9f9; }
.summary-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 15px; }
.summary-card { background: #ecf0f1; padding: 15px; border-radius: 8px; text-align: center; }
.summary-card .value { font-size: 28px; font-weight: bold; color: #2c3e50; }
.summary-card .label { font-size: 13px; color: #7f8c8d; }
.controls { display: flex; gap: 15px; align-items: center; margin-bottom: 10px; flex-wrap: wrap; }
.search-box { padding: 6px 12px; border: 1px solid #bdc3c7; border-radius: 4px; font-size: 13px; width: 250px; }
.page-size-select { padding: 6px; border: 1px solid #bdc3c7; border-radius: 4px; font-size: 13px; }
.btn { padding: 8px 16px; border: none; border-radius: 4px; cursor: pointer; font-size: 13px; }
.btn-primary { background-color: #3498db; color: white; }
.btn-primary:hover { background-color: #2980b9; }
.btn-secondary { background-color: #ecf0f1; color: #2c3e50; }
.btn-secondary:hover { background-color: #bdc3c7; }
.pagination { margin-top: 10px; display: flex; gap: 5px; align-items: center; }
.pagination span { font-size: 13px; color: #7f8c8d; }
.seq-cell { font-family: monospace; font-size: 11px; max-width: 300px; overflow: hidden; text-overflow: ellipsis; }
.footer { margin-top: 40px; padding-top: 15px; border-top: 1px solid #ecf0f1; font-size: 12px; color: #95a5a6; }
.legend { display: flex; gap: 15px; flex-wrap: wrap; margin-bottom: 10px; font-size: 12px; }
.legend-item { display: flex; align-items: center; gap: 5px; }
.legend-dot { width: 12px; height: 12px; border-radius: 50%; display: inline-block; }
.annotated-seq { background: #fafafa; padding: 15px; border: 1px solid #ddd; border-radius: 4px; overflow-x: auto; }
.seq-legend { display: flex; gap: 20px; flex-wrap: wrap; margin-bottom: 10px; font-size: 12px; }
.seq-legend-item { display: flex; align-items: center; gap: 5px; }
.seq-legend-box { width: 20px; height: 14px; display: inline-block; border-radius: 2px; }
.mutation-tip { cursor: help; position: relative; }
.mutation-tip:hover::after {
  content: attr(data-tip);
  position: absolute; bottom: 100%; left: 50%; transform: translateX(-50%);
  background: #2c3e50; color: white; padding: 6px 10px; border-radius: 4px;
  font-size: 11px; white-space: nowrap; z-index: 100; pointer-events: none;
  font-weight: normal; font-family: 'Segoe UI', Arial, sans-serif;
  margin-bottom: 4px; max-width: 400px; white-space: normal; line-height: 1.4;
}
.s { font-family: 'Courier New', monospace; font-size: 12px; }
.su { background: #e0e0e0; }
.sm { font-weight: bold; }
.sm-p { color: #e74c3c; }
.sm-lp { color: #e67e22; }
.sm-v { color: #f39c12; }
.sm-b { color: #2ecc71; }
.sm-lb { color: #27ae60; }
.sm-o { color: #95a5a6; }
.sg { font-family: 'Courier New', monospace; }
.ll { font-family: 'Courier New', monospace; font-size: 12px; color: #2980b9; }
.lp { font-family: 'Courier New', monospace; font-size: 12px; color: transparent; user-select: none; }
.sl { height: 14px; white-space: pre; line-height: 1; }
.sn { white-space: pre; }
.se { margin: 15px 0; }
.sh { font-weight: bold; color: #2c3e50; margin-bottom: 4px; }
.sp { height: 14px; }
.summary-info { margin-top: 15px; color: #7f8c8d; }
.table-wrap { overflow-x: auto; }
.copy-hidden { display: none; }
.ld-p { background: #e74c3c; }
.ld-lp { background: #e67e22; }
.ld-v { background: #f39c12; }
.ld-b { background: #2ecc71; }
.ld-pr { background: rgba(52,152,219,0.5); }
.slb-u { background: #e0e0e0; }
.slb-po { background: #dee6ef; }
.slb-pe { background: #ffffd7; }
.slb-ov { background: #eef2e3; }
.slm { color: #e74c3c; font-weight: bold; font-family: monospace; }
</style>
</head>
<body>

<h1>{{ gene }} Probe Design Report</h1>

<div class="section">
  <h2>Summary</h2>
  <div class="summary-grid">
    <div class="summary-card">
      <div class="value">{{ n_mutations }}</div>
      <div class="label">ClinVar Missense Mutations</div>
    </div>
    <div class="summary-card">
      <div class="value">{{ n_probes }}</div>
      <div class="label">Probes Designed</div>
    </div>
    <div class="summary-card">
      <div class="value">{{ coverage_pct }}%</div>
      <div class="label">Mutation Coverage ({{ n_covered }}/{{ n_unique_pos }})</div>
    </div>
  </div>
  <p class="summary-info">
    Transcript: <b>{{ accession }}</b> &nbsp;|&nbsp;
    Sequence length: <b>{{ seq_len }}</b> nt &nbsp;|&nbsp;
    CDS: <b>{{ cds_start }}-{{ cds_end }}</b> &nbsp;|&nbsp;
    Exons: <b>{{ n_exons }}</b> &nbsp;|&nbsp;
    Junctions: <b>{{ n_junctions }}</b>
  </p>
</div>

<div class="section">
  <h2>Transcript Map</h2>
  <div class="legend">
    <div class="legend-item"><span class="legend-dot ld-p"></span> Pathogenic</div>
    <div class="legend-item"><span class="legend-dot ld-lp"></span> Likely pathogenic</div>
    <div class="legend-item"><span class="legend-dot ld-v"></span> VUS</div>
    <div class="legend-item"><span class="legend-dot ld-b"></span> Benign/Likely benign</div>
    <div class="legend-item"><span class="legend-dot ld-pr"></span> Probes</div>
  </div>
  <div id="plotly-chart"></div>
  <script>
    var figure = {{ plotly_json|safe }};
    Plotly.newPlot('plotly-chart', figure.data, figure.layout);
  </script>
</div>

<div class="section">
  <h2>ClinVar Mutations</h2>
  <div class="controls">
    <input type="text" id="mutation-search" placeholder="Search mutations..." class="search-box">
    <label>Show
      <select id="mutation-page-size" class="page-size-select">
        <option value="10" selected>10</option>
        <option value="25">25</option>
        <option value="50">50</option>
        <option value="100">100</option>
        <option value="-1">All</option>
      </select>
      entries
    </label>
    <button class="btn btn-primary" onclick="downloadMutationsCSV()">Download CSV</button>
  </div>
  <div class="table-wrap">
  <table id="mutation-table">
    <thead>
      <tr>
        <th data-col="id">ID</th>
        <th data-col="protein_change">Protein Change</th>
        <th data-col="cdna_change">cDNA Change</th>
        <th data-col="clinical_significance">Significance</th>
        <th data-col="mrna_position">mRNA Position</th>
        <th data-col="cdna_position">CDS Position</th>
        <th data-col="clinvar_accession">ClinVar Accession</th>
        <th data-col="chrom">Chrom</th>
        <th data-col="pos_start">Genomic Pos</th>
        <th data-col="exon">Exon</th>
        <th data-col="probe">Probe</th>
      </tr>
    </thead>
    <tbody id="mutation-tbody"></tbody>
  </table>
  </div>
  <div class="pagination" id="mutation-pagination"></div>
</div>

<div class="section">
  <h2>Designed Probes</h2>
  <div class="controls">
    <input type="text" id="probe-search" placeholder="Search probes..." class="search-box">
    <label>Show
      <select id="probe-page-size" class="page-size-select">
        <option value="10" selected>10</option>
        <option value="25">25</option>
        <option value="50">50</option>
        <option value="100">100</option>
        <option value="-1">All</option>
      </select>
      entries
    </label>
    <button class="btn btn-primary" onclick="downloadProbesCSV()">Download CSV</button>
  </div>
  <div class="table-wrap">
  <table id="probe-table">
    <thead>
      <tr>
        <th data-col="probe_name">Probe Name</th>
        <th data-col="exon">Exon</th>
        <th data-col="mrna_start">mRNA Start</th>
        <th data-col="mrna_end">mRNA End</th>
        <th data-col="cds_start">CDS Start</th>
        <th data-col="cds_end">CDS End</th>
        <th data-col="length">Length</th>
        <th data-col="target_sequence">Target</th>
        <th data-col="fwd_oligo">Forward Oligo (5'->3')</th>
        <th data-col="rev_oligo">Reverse Oligo (5'->3')</th>
        <th data-col="mutations_covered">Mutations Covered</th>
        <th data-col="cdna_changes_covered">cDNA Changes</th>
      </tr>
    </thead>
    <tbody id="probe-tbody"></tbody>
  </table>
  </div>
  <div class="pagination" id="probe-pagination"></div>
</div>

<div class="section">
  <h2>Annotated mRNA Sequence</h2>
  <div class="seq-legend">
    <div class="seq-legend-item"><span class="seq-legend-box slb-u"></span> 5'/3' UTR</div>
    <div class="seq-legend-item"><span class="slm">N</span> Mutation</div>
    <div class="seq-legend-item"><span class="seq-legend-box slb-po"></span> Probe (odd)</div>
    <div class="seq-legend-item"><span class="seq-legend-box slb-pe"></span> Probe (even)</div>
    <div class="seq-legend-item"><span class="seq-legend-box slb-ov"></span> Probe overlap</div>
  </div>
  <div class="controls">
    <button class="btn btn-primary" onclick="copyAnnotatedSequence()">Copy Sequence</button>
  </div>
  <div id="annotated-sequence" class="annotated-seq">
    {{ annotated_html|safe }}
  </div>
  <div id="copy-sequence" class="copy-hidden">
    {{ copy_html|safe }}
  </div>
</div>

<div class="footer">Generated by gene-probe-design pipeline on {{ date }}</div>

<script>
const mutationsData = {{ mutations_json|safe }};
const probesData = {{ probes_json|safe }};

let mutPage = 0;
let mutPageSize = 10;
let mutSortCol = null;
let mutSortAsc = true;
let filteredMutations = [...mutationsData];

let probePage = 0;
let probePageSize = 10;
let probeSortCol = null;
let probeSortAsc = true;
let filteredProbes = [...probesData];

function renderMutations() {
  const tbody = document.getElementById('mutation-tbody');
  const start = mutPageSize === -1 ? 0 : mutPage * mutPageSize;
  const end = mutPageSize === -1 ? filteredMutations.length : start + mutPageSize;
  const page = filteredMutations.slice(start, end);

  tbody.innerHTML = page.map(m => `<tr>
    <td>${m.id || ''}</td>
    <td>${m.protein_change || ''}</td>
    <td>${m.cdna_change || ''}</td>
    <td>${m.clinical_significance || ''}</td>
    <td>${m.mrna_position != null ? m.mrna_position : ''}</td>
    <td>${m.cdna_position != null ? m.cdna_position : ''}</td>
    <td>${m.clinvar_accession || ''}</td>
    <td>${m.chrom || ''}</td>
    <td>${m.pos_start || ''}</td>
    <td>${m.exon || ''}</td>
    <td>${m.probe || ''}</td>
  </tr>`).join('');

  const totalPages = mutPageSize === -1 ? 1 : Math.ceil(filteredMutations.length / mutPageSize);
  const pag = document.getElementById('mutation-pagination');
  pag.innerHTML = `<span>Showing ${start+1}-${Math.min(end, filteredMutations.length)} of ${filteredMutations.length}</span> ` +
    (mutPage > 0 ? `<button class="btn btn-secondary" onclick="mutPage--;renderMutations()">Prev</button>` : '') +
    (mutPage < totalPages - 1 ? `<button class="btn btn-secondary" onclick="mutPage++;renderMutations()">Next</button>` : '');
}

function renderProbes() {
  const tbody = document.getElementById('probe-tbody');
  const start = probePageSize === -1 ? 0 : probePage * probePageSize;
  const end = probePageSize === -1 ? filteredProbes.length : start + probePageSize;
  const page = filteredProbes.slice(start, end);

  tbody.innerHTML = page.map(p => `<tr>
    <td>${p.probe_name}</td>
    <td>${p.exon || ''}</td>
    <td>${p.mrna_start}</td><td>${p.mrna_end}</td>
    <td>${p.cds_start}</td><td>${p.cds_end}</td>
    <td>${p.length}</td>
    <td class="seq-cell">${p.target_sequence}</td>
    <td class="seq-cell">${p.fwd_oligo}</td>
    <td class="seq-cell">${p.rev_oligo}</td>
    <td>${p.mutations_covered}</td>
    <td>${p.cdna_changes_covered}</td>
  </tr>`).join('');

  const totalPages = probePageSize === -1 ? 1 : Math.ceil(filteredProbes.length / probePageSize);
  const pag = document.getElementById('probe-pagination');
  pag.innerHTML = `<span>Showing ${start+1}-${Math.min(end, filteredProbes.length)} of ${filteredProbes.length}</span> ` +
    (probePage > 0 ? `<button class="btn btn-secondary" onclick="probePage--;renderProbes()">Prev</button>` : '') +
    (probePage < totalPages - 1 ? `<button class="btn btn-secondary" onclick="probePage++;renderProbes()">Next</button>` : '');
}

function filterProbes() {
  const q = document.getElementById('probe-search').value.toLowerCase();
  filteredProbes = probesData.filter(p =>
    Object.values(p).some(v => String(v).toLowerCase().includes(q))
  );
  probePage = 0;
  renderProbes();
}

function sortProbes(col) {
  if (probeSortCol === col) { probeSortAsc = !probeSortAsc; }
  else { probeSortCol = col; probeSortAsc = true; }
  filteredProbes.sort((a, b) => {
    let va = a[col], vb = b[col];
    if (va == null) return 1;
    if (vb == null) return -1;
    if (typeof va === 'number') return probeSortAsc ? va - vb : vb - va;
    return probeSortAsc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
  });
  probePage = 0;
  renderProbes();
}

function filterMutations() {
  const q = document.getElementById('mutation-search').value.toLowerCase();
  filteredMutations = mutationsData.filter(m =>
    Object.values(m).some(v => String(v).toLowerCase().includes(q))
  );
  mutPage = 0;
  renderMutations();
}

function sortMutations(col) {
  if (mutSortCol === col) { mutSortAsc = !mutSortAsc; }
  else { mutSortCol = col; mutSortAsc = true; }
  filteredMutations.sort((a, b) => {
    let va = a[col], vb = b[col];
    if (va == null) return 1;
    if (vb == null) return -1;
    if (typeof va === 'number') return mutSortAsc ? va - vb : vb - va;
    return mutSortAsc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
  });
  mutPage = 0;
  renderMutations();
}

function downloadMutationsCSV() {
  const header = 'ID,Protein Change,cDNA Change,Significance,mRNA Position,CDS Position,ClinVar Accession,Chrom,Genomic Pos,Exon,Probe';
  const rows = filteredMutations.map(m =>
    [m.id||'', m.protein_change||'', m.cdna_change||'', m.clinical_significance||'',
     m.mrna_position!=null?m.mrna_position:'', m.cdna_position!=null?m.cdna_position:'',
     m.clinvar_accession||'', m.chrom||'', m.pos_start||'', m.exon||'', m.probe||'']
    .map(v => '"'+String(v).replace(/"/g,'""')+'"').join(',')
  );
  const csv = header + '\\n' + rows.join('\\n');
  const blob = new Blob([csv], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = '{{ gene }}_mutations.csv';
  a.click();
}

function downloadProbesCSV() {
  const header = 'probe_name,exon,FWD_name,REV_name,mRNA_start,mRNA_end,CDS_start,CDS_end,length,target_sequence,fwd_oligo,rev_oligo,mutations_covered,cdna_changes_covered,mutation_positions';
  const rows = probesData.map(p =>
    [p.probe_name, p.exon||'', p.probe_name+'-FWD', p.probe_name+'-REV',
     p.mrna_start, p.mrna_end, p.cds_start, p.cds_end, p.length,
     p.target_sequence, p.fwd_oligo, p.rev_oligo,
     p.mutations_covered, p.cdna_changes_covered, p.mutation_positions]
    .map(v => '"'+String(v).replace(/"/g,'""')+'"').join(',')
  );
  const csv = header + '\\n' + rows.join('\\n');
  const blob = new Blob([csv], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = '{{ gene }}_probes.csv';
  a.click();
}

function copyAnnotatedSequence() {
  const el = document.getElementById('copy-sequence');
  // Temporarily show the hidden div so innerText works correctly
  el.style.display = 'block';
  el.style.position = 'absolute';
  el.style.left = '-9999px';
  const html = el.innerHTML;
  const text = el.innerText;
  el.style.display = 'none';
  el.style.position = '';
  el.style.left = '';
  const blob = new Blob([html], { type: 'text/html' });
  const textBlob = new Blob([text], { type: 'text/plain' });
  if (navigator.clipboard && navigator.clipboard.write) {
    navigator.clipboard.write([
      new ClipboardItem({ 'text/html': blob, 'text/plain': textBlob })
    ]).then(() => alert('Sequence copied with formatting!')).catch(() => {
      const range = document.createRange();
      range.selectNodeContents(el);
      window.getSelection().removeAllRanges();
      window.getSelection().addRange(range);
      document.execCommand('copy');
      alert('Sequence copied!');
    });
  } else {
    el.style.display = 'block';
    el.style.position = 'absolute';
    el.style.left = '-9999px';
    const range = document.createRange();
    range.selectNodeContents(el);
    window.getSelection().removeAllRanges();
    window.getSelection().addRange(range);
    document.execCommand('copy');
    el.style.display = 'none';
    el.style.position = '';
    el.style.left = '';
    alert('Sequence copied!');
  }
}

document.getElementById('mutation-search').addEventListener('input', filterMutations);
document.getElementById('mutation-page-size').addEventListener('change', function() {
  mutPageSize = parseInt(this.value);
  mutPage = 0;
  renderMutations();
});
document.querySelectorAll('#mutation-table th').forEach(th => {
  th.addEventListener('click', () => sortMutations(th.dataset.col));
});

document.getElementById('probe-search').addEventListener('input', filterProbes);
document.getElementById('probe-page-size').addEventListener('change', function() {
  probePageSize = parseInt(this.value);
  probePage = 0;
  renderProbes();
});
document.querySelectorAll('#probe-table th').forEach(th => {
  th.addEventListener('click', () => sortProbes(th.dataset.col));
});

renderMutations();
renderProbes();
</script>
</body>
</html>
"""


def main():
    gene = snakemake.params.gene
    region_config = snakemake.params.region

    mutations_df, sequence, accession, cds_start, cds_end, junctions_df, probes_df = load_data(
        snakemake.input.mutations,
        snakemake.input.transcript,
        snakemake.input.junctions,
        snakemake.input.probes,
    )

    seq_len = len(sequence)
    n_exons = len(junctions_df) + 1

    # Enrich mutations with computed fields
    enriched_mutations = build_enriched_mutations(
        mutations_df, cds_start, cds_end, junctions_df, probes_df,
    )

    # Filter mutations by region config
    if region_config and region_config != "All":
        junction_positions = junctions_df["mrna_position"].tolist()
        region_tokens = _parse_region_config(region_config)
        allowed_ranges = _mrna_ranges_for_regions(
            region_tokens, junction_positions, seq_len, cds_start, cds_end,
        )

        def _in_region(mrna_pos):
            if pd.isna(mrna_pos):
                return False
            return any(s <= mrna_pos <= e for s, e in allowed_ranges)

        n_before = len(enriched_mutations)
        enriched_mutations = enriched_mutations[
            enriched_mutations["mrna_position"].apply(_in_region)
        ]
        n_after = len(enriched_mutations)
        if n_after < n_before:
            log.info(
                f"Region filter '{region_config}': kept {n_after}/{n_before} mutations in report"
            )

    # Compute coverage stats (only mutations with valid mRNA positions)
    valid_mut = enriched_mutations[enriched_mutations["mrna_position"].notna()]
    all_positions = set(valid_mut["mrna_position"].astype(int).tolist())
    covered_positions = set()
    for _, row in probes_df.iterrows():
        for p in str(row["mutation_positions"]).split(","):
            if p.strip():
                covered_positions.add(int(p.strip()))
    n_unique_pos = len(all_positions)
    n_covered = len(covered_positions & all_positions)
    coverage_pct = f"{100 * n_covered / n_unique_pos:.1f}" if n_unique_pos else "0"

    log.info(
        f"Report for {gene}: {len(enriched_mutations)} mutations, {len(probes_df)} probes, "
        f"{coverage_pct}% coverage"
    )

    # Build Plotly figure
    fig = build_plotly_figure(enriched_mutations, accession, junctions_df, probes_df, seq_len)
    plotly_json = fig.to_json()

    # Build annotated sequence HTML
    annotated_html, copy_html = build_annotated_sequence(
        sequence, cds_start, cds_end, junctions_df, enriched_mutations, probes_df
    )

    # Serialize data for JavaScript
    mutations_json = enriched_mutations.to_json(orient="records")
    probes_json = probes_df.to_json(orient="records")

    # Render HTML
    html = Template(HTML_TEMPLATE).render(
        gene=gene,
        accession=accession,
        seq_len=seq_len,
        cds_start=cds_start,
        cds_end=cds_end,
        n_mutations=len(enriched_mutations),
        n_probes=len(probes_df),
        n_exons=n_exons,
        n_junctions=len(junctions_df),
        n_covered=n_covered,
        n_unique_pos=n_unique_pos,
        coverage_pct=coverage_pct,
        plotly_json=plotly_json,
        mutations_json=mutations_json,
        probes_json=probes_json,
        annotated_html=annotated_html,
        copy_html=copy_html,
        date=date.today().isoformat(),
    )

    with open(snakemake.output.report, "w") as f:
        f.write(html)
    log.info(f"Wrote report to {snakemake.output.report}")

    # Write detailed mutations file
    write_mutations_file(enriched_mutations, snakemake.output.mutations_file)
    log.info(f"Wrote mutations file to {snakemake.output.mutations_file}")

    # Write IDT order file
    if hasattr(snakemake.output, "idt_order") and snakemake.output.idt_order:
        idt_rows = []
        for _, prow in probes_df.iterrows():
            idt_rows.append({
                "Probe_name_fwd": f"{prow['probe_name']}-FWD",
                "fwd_seq": prow["fwd_oligo"],
                "Probe_name_rev": f"{prow['probe_name']}-REV",
                "rev_seq": prow["rev_oligo"],
            })
        idt_df = pd.DataFrame(idt_rows)
        idt_df.to_csv(snakemake.output.idt_order, sep="\t", index=False)
        log.info(f"Wrote IDT order file to {snakemake.output.idt_order}")


if __name__ == "__main__":
    main()
