"""Generate interactive HTML report with Plotly visualization."""

import json
import logging

import pandas as pd
import plotly.graph_objects as go
from jinja2 import Template

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def load_data(mutations_path, transcript_path, junctions_path, probes_path):
    mutations_df = pd.read_csv(mutations_path, sep="\t")
    transcript_df = pd.read_csv(transcript_path, sep="\t")
    junctions_df = pd.read_csv(junctions_path, sep="\t")
    probes_df = pd.read_csv(probes_path, sep="\t")

    sequence = transcript_df.loc[
        transcript_df["field"] == "sequence", "value"
    ].values[0]
    accession = transcript_df.loc[
        transcript_df["field"] == "accession", "value"
    ].values[0]

    return mutations_df, sequence, accession, junctions_df, probes_df


def build_plotly_figure(
    mutations_df, accession, junctions_df, probes_df, seq_len
):
    fig = go.Figure()

    # --- Exon blocks ---
    exons = []
    prev_end = 0
    for _, jrow in junctions_df.iterrows():
        exons.append((prev_end + 1, jrow["mrna_position"]))
        prev_end = jrow["mrna_position"]
    exons.append((prev_end + 1, seq_len))

    colors = []
    n_exons = len(exons)
    for i in range(n_exons):
        colors.append(f"hsl({int(210 + 40 * (i / n_exons))}, 60%, 75%)")

    for i, (start, end) in enumerate(exons):
        fig.add_trace(
            go.Bar(
                x=[end - start + 1],
                y=[1],
                base=start - 1,
                orientation="h",
                marker_color=colors[i],
                name=f"Exon {i + 1}",
                hovertemplate=f"Exon {i + 1}<br>Position: {start}-{end}<extra></extra>",
                showlegend=False,
            )
        )

    # --- Mutations ---
    sig_colors = {
        "Pathogenic": "#e74c3c",
        "Pathogenic/Likely pathogenic": "#e74c3c",
        "Likely pathogenic": "#e67e22",
        "Uncertain significance": "#f39c12",
        "Likely benign": "#27ae60",
        "Benign": "#2ecc71",
        "Benign/Likely benign": "#2ecc71",
    }

    mut_colors = [
        sig_colors.get(s, "#95a5a6") for s in mutations_df["clinical_significance"]
    ]

    hover_text = [
        f"<b>{row['protein_change']}</b><br>"
        f"cDNA: {row['cdna_change']}<br>"
        f"Significance: {row['clinical_significance']}<br>"
        f"Accession: {row['clinvar_accession']}<br>"
        f"Position: {row['mrna_position']}"
        for _, row in mutations_df.iterrows()
    ]

    fig.add_trace(
        go.Scatter(
            x=mutations_df["mrna_position"].tolist(),
            y=[1.8] * len(mutations_df),
            mode="markers",
            marker=dict(size=6, color=mut_colors, line=dict(width=0.5, color="white")),
            text=hover_text,
            hoverinfo="text",
            name="Mutations",
            showlegend=False,
        )
    )

    # --- Probes ---
    for _, prow in probes_df.iterrows():
        muts = prow["mutations_covered"].replace(";", ", ")
        fig.add_trace(
            go.Bar(
                x=[prow["end"] - prow["start"] + 1],
                y=[2.6],
                base=prow["start"] - 1,
                orientation="h",
                marker_color="rgba(52, 152, 219, 0.5)",
                marker_line=dict(color="rgba(41, 128, 185, 0.8)", width=1),
                name=prow["probe_name"],
                hovertemplate=(
                    f"<b>{prow['probe_name']}</b><br>"
                    f"Position: {prow['start']}-{prow['end']} ({prow['length']} nt)<br>"
                    f"Mutations: {muts}<extra></extra>"
                ),
                showlegend=False,
            )
        )

    fig.update_layout(
        barmode="overlay",
        xaxis=dict(title="mRNA Position (nt)", rangeslider=dict(visible=True)),
        yaxis=dict(
            tickvals=[1, 1.8, 2.6],
            ticktext=["Exons", "Mutations", "Probes"],
            range=[0.3, 3.2],
        ),
        height=400,
        margin=dict(l=80, r=30, t=30, b=60),
        plot_bgcolor="white",
    )

    return fig


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
  <p style="margin-top:15px; color:#7f8c8d;">
    Transcript: <b>{{ accession }}</b> &nbsp;|&nbsp;
    Sequence length: <b>{{ seq_len }}</b> nt &nbsp;|&nbsp;
    Exons: <b>{{ n_exons }}</b> &nbsp;|&nbsp;
    Junctions: <b>{{ n_junctions }}</b>
  </p>
</div>

<div class="section">
  <h2>Transcript Map</h2>
  <div class="legend">
    <div class="legend-item"><span class="legend-dot" style="background:#e74c3c;"></span> Pathogenic</div>
    <div class="legend-item"><span class="legend-dot" style="background:#e67e22;"></span> Likely pathogenic</div>
    <div class="legend-item"><span class="legend-dot" style="background:#f39c12;"></span> VUS</div>
    <div class="legend-item"><span class="legend-dot" style="background:#2ecc71;"></span> Benign/Likely benign</div>
    <div class="legend-item"><span class="legend-dot" style="background:rgba(52,152,219,0.5);"></span> Probes</div>
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
  </div>
  <table id="mutation-table">
    <thead>
      <tr>
        <th data-col="protein_change">Protein Change</th>
        <th data-col="cdna_change">cDNA Change</th>
        <th data-col="clinical_significance">Significance</th>
        <th data-col="mrna_position">mRNA Position</th>
        <th data-col="clinvar_accession">ClinVar Accession</th>
        <th data-col="chrom">Chrom</th>
        <th data-col="pos_start">Genomic Pos</th>
      </tr>
    </thead>
    <tbody id="mutation-tbody"></tbody>
  </table>
  <div class="pagination" id="mutation-pagination"></div>
</div>

<div class="section">
  <h2>Designed Probes</h2>
  <div class="controls">
    <button class="btn btn-primary" onclick="downloadCSV()">Download CSV</button>
  </div>
  <table id="probe-table">
    <thead>
      <tr>
        <th>Probe Name</th>
        <th>Start</th>
        <th>End</th>
        <th>Length</th>
        <th>Forward Oligo (5'->3')</th>
        <th>Reverse Oligo (5'->3')</th>
        <th>Mutations Covered</th>
      </tr>
    </thead>
    <tbody id="probe-tbody"></tbody>
  </table>
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

function renderMutations() {
  const tbody = document.getElementById('mutation-tbody');
  const start = mutPageSize === -1 ? 0 : mutPage * mutPageSize;
  const end = mutPageSize === -1 ? filteredMutations.length : start + mutPageSize;
  const page = filteredMutations.slice(start, end);

  tbody.innerHTML = page.map(m => `<tr>
    <td>${m.protein_change}</td><td>${m.cdna_change}</td>
    <td>${m.clinical_significance}</td><td>${m.mrna_position}</td>
    <td>${m.clinvar_accession}</td><td>${m.chrom}</td><td>${m.pos_start}</td>
  </tr>`).join('');

  const totalPages = mutPageSize === -1 ? 1 : Math.ceil(filteredMutations.length / mutPageSize);
  const pag = document.getElementById('mutation-pagination');
  pag.innerHTML = `<span>Showing ${start+1}-${Math.min(end, filteredMutations.length)} of ${filteredMutations.length}</span> ` +
    (mutPage > 0 ? `<button class="btn btn-secondary" onclick="mutPage--;renderMutations()">Prev</button>` : '') +
    (mutPage < totalPages - 1 ? `<button class="btn btn-secondary" onclick="mutPage++;renderMutations()">Next</button>` : '');
}

function renderProbes() {
  document.getElementById('probe-tbody').innerHTML = probesData.map(p => `<tr>
    <td>${p.probe_name}</td><td>${p.start}</td><td>${p.end}</td><td>${p.length}</td>
    <td class="seq-cell">${p.fwd_oligo}</td><td class="seq-cell">${p.rev_oligo}</td>
    <td>${p.mutations_covered.replace(/;/g, ', ')}</td>
  </tr>`).join('');
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
    if (typeof va === 'number') return mutSortAsc ? va - vb : vb - va;
    return mutSortAsc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
  });
  mutPage = 0;
  renderMutations();
}

function downloadCSV() {
  const header = 'probe_name,FWD_name,REV_name,start,end,length,fwd_oligo,rev_oligo,mutations_covered,mutation_positions';
  const rows = probesData.map(p =>
    `${p.probe_name},${p.probe_name}-FWD,${p.probe_name}-REV,${p.start},${p.end},${p.length},${p.fwd_oligo},${p.rev_oligo},${p.mutations_covered},${p.mutation_positions}`
  );
  const csv = header + '\\n' + rows.join('\\n');
  const blob = new Blob([csv], { type: 'text/csv' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = '{{ gene }}_probes.csv';
  a.click();
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

renderMutations();
renderProbes();
</script>
</body>
</html>
"""


def main():
    gene = snakemake.params.gene

    mutations_df, sequence, accession, junctions_df, probes_df = load_data(
        snakemake.input.mutations,
        snakemake.input.transcript,
        snakemake.input.junctions,
        snakemake.input.probes,
    )

    seq_len = len(sequence)
    n_exons = len(junctions_df) + 1

    # Compute coverage stats
    all_positions = set(mutations_df["mrna_position"].tolist())
    covered_positions = set()
    for _, row in probes_df.iterrows():
        for p in str(row["mutation_positions"]).split(";"):
            covered_positions.add(int(p))
    n_unique_pos = len(all_positions)
    n_covered = len(covered_positions)
    coverage_pct = f"{100 * n_covered / n_unique_pos:.1f}" if n_unique_pos else "0"

    log.info(
        f"Report for {gene}: {len(mutations_df)} mutations, {len(probes_df)} probes, "
        f"{coverage_pct}% coverage"
    )

    fig = build_plotly_figure(
        mutations_df, accession, junctions_df, probes_df, seq_len
    )
    plotly_json = fig.to_json()

    mutations_json = mutations_df.to_json(orient="records")
    probes_json = probes_df.to_json(orient="records")

    from datetime import date
    html = Template(HTML_TEMPLATE).render(
        gene=gene,
        accession=accession,
        seq_len=seq_len,
        n_mutations=len(mutations_df),
        n_probes=len(probes_df),
        n_exons=n_exons,
        n_junctions=len(junctions_df),
        n_covered=n_covered,
        n_unique_pos=n_unique_pos,
        coverage_pct=coverage_pct,
        plotly_json=plotly_json,
        mutations_json=mutations_json,
        probes_json=probes_json,
        date=date.today().isoformat(),
    )

    with open(snakemake.output[0], "w") as f:
        f.write(html)

    log.info(f"Wrote report to {snakemake.output[0]}")


if __name__ == "__main__":
    main()
