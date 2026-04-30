# Gene Probe Design Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A57.0-green.svg)](https://snakemake.readthedocs.io/)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)

A Snakemake pipeline for designing antisense gene probes and oligos. Given a gene name in `config.yaml`, the pipeline automatically fetches ClinVar variants (configurable molecular consequence types) and transcript data from NCBI, designs probes (24--26 nt) that maximize mutation coverage while avoiding exon-exon junction buffer zones, and produces an interactive HTML report with orderable oligo sequences.

Designed for any gene --- POLR3A is the example, but changing the gene list in `config.yaml` is all that's needed to target a different gene.

**Key features:**

- Works with any gene --- just change `config.yaml`
- Fully programmatic data fetching from NCBI (ClinVar + GenBank)
- Greedy probe placement algorithm with junction-aware avoidance
- Interactive HTML report with Plotly visualizations, sortable tables, and annotated sequence
- Order-ready oligo sequences (forward and reverse)

## Pipeline Overview

```
config.yaml
    │
    ▼
[fetch_clinvar] ──► data/{gene}/clinvar_mutations.tsv
    │                 data/{gene}/principal_transcript.txt
    ▼
[fetch_transcript] ──► data/{gene}/transcript.tsv
    │                    data/{gene}/junctions.tsv
    ▼
[design_probes] ──► data/{gene}/probes.tsv
    │
    ├──► [generate_report] ──► results/{gene}/{gene}_report.html
    │                          results/{gene}/{gene}_mutations.tsv
    └──► [export_probes]     ──► results/{gene}/probes.tsv
```

## Installation

### Prerequisites

- **Conda**: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- **Snakemake** >= 7.0 (installed via conda, see below)

### Linux / macOS

1. **Install Miniconda** (if you don't have conda already):

   Download and run the installer from [docs.conda.io](https://docs.conda.io/en/latest/miniconda.html).

2. **Create and activate a Snakemake environment:**

   ```bash
   conda create -n snakemake7.24 -c conda-forge -c bioconda snakemake=7.24
   conda activate snakemake7.24
   ```

3. **Clone this repository:**

   ```bash
   git clone <repository-url>
   cd pipeline
   ```

4. **Verify the installation with a dry run:**

   ```bash
   snakemake -n --cores 2
   ```

### Windows

Snakemake does not have official Windows support. There are two options:

#### Option A: WSL2 (recommended)

1. Install [Windows Subsystem for Linux (WSL2)](https://learn.microsoft.com/en-us/windows/wsl/install) with Ubuntu:

   ```powershell
   wsl --install
   ```

2. Open a WSL terminal and follow the **Linux / macOS** instructions above.

#### Option B: Native Windows (experimental)

1. Install [Miniconda for Windows](https://docs.conda.io/en/latest/miniconda.html).
2. Open Anaconda Prompt and create the environment:

   ```cmd
   conda create -n snakemake7.24 -c conda-forge -c bioconda snakemake=7.24
   conda activate snakemake7.24
   ```

3. Clone and run --- note that some Snakemake features may not work correctly on Windows. WSL2 is strongly recommended.

## Quick Start

1. **Edit `config.yaml`** with your gene of interest:

   ```yaml
   genes:
     - POLR3A
   ```

2. **Run the pipeline:**

   ```bash
   conda activate snakemake7.24
   snakemake --use-conda --cores 2
   ```

   On the first run, Snakemake will create conda environments for each rule automatically (this may take a few minutes). Subsequent runs reuse them.

3. **View the results:**

   ```bash
   # Open the interactive HTML report in your browser
   xdg-open results/POLR3A/POLR3A_report.html     # Linux
   open results/POLR3A/POLR3A_report.html           # macOS
   ```

   Output files are in `results/{gene}/`:
   - `POLR3A_report.html` --- interactive report with transcript map, tables, and annotated sequence
   - `probes.tsv` --- probe sequences and oligo designs for ordering
   - `POLR3A_mutations.tsv` --- detailed mutation data

## Configuration

All settings are in `config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `genes` | `[POLR3A]` | List of gene symbols. The pipeline runs independently for each gene. |
| `ncbi_api_key` | `""` | Optional NCBI API key. Increases rate limit from 3 to 10 requests/s. [Get one here](https://www.ncbi.nlm.nih.gov/account/settings/). |
| `molecular_consequences` | `["missense variant"]` | ClinVar molecular consequence filter. Accepts a list of values (e.g. `["missense variant", "nonsense"]`), a single value, or `"all"` to fetch every type. See below for valid values. |
| `probe_params.target_length` | `25` | Ideal probe size in nucleotides |
| `probe_params.min_length` | `24` | Minimum acceptable probe size |
| `probe_params.max_length` | `26` | Maximum acceptable probe size |
| `probe_params.junction_buffer` | `5` | Nucleotides excluded on each side of exon-exon junctions |
| `probe_params.scan_offset` | `5` | Start scanning this many nt from the 5' end |

**To analyze a different gene**, simply change the gene list:

```yaml
genes:
  - BRCA1
  - TP53
```

The principal transcript is auto-discovered from ClinVar data (most common NM\_ accession across all fetched variant records). No manual transcript ID is needed.

### Molecular Consequence Values

These are the valid values for `molecular_consequences` in `config.yaml`:

| Value | Description |
|-------|-------------|
| `missense variant` | Amino acid substitution (default) |
| `nonsense` | Stop codon gained |
| `frameshift variant` | Insertion/deletion causing frameshift |
| `synonymous variant` | No amino acid change |
| `intron variant` | Within intron |
| `non-coding transcript variant` | In non-coding transcript |
| `5 prime UTR variant` | In 5' UTR |
| `3 prime UTR variant` | In 3' UTR |
| `splice donor variant` | Splice donor site |
| `splice acceptor variant` | Splice acceptor site |
| `inframe deletion` | In-frame deletion |
| `inframe insertion` | In-frame insertion |
| `inframe indel` | In-frame insertion/deletion |
| `initiator codon variant` | Start codon change |
| `stop lost` | Stop codon lost |
| `genic upstream transcript variant` | Upstream transcript overlap |
| `genic downstream transcript variant` | Downstream transcript overlap |
| `no sequence alteration` | No sequence change |

Use `"all"` to fetch variants of every consequence type (no filter applied).

**Example:** Fetch both missense and nonsense variants:

```yaml
molecular_consequences:
  - missense variant
  - nonsense
```

## Output Files

For each gene, the pipeline produces three output files in `results/{gene}/`:

| File | Description |
|------|-------------|
| `{gene}_report.html` | Self-contained interactive HTML report with Plotly visualizations, sortable/searchable mutation and probe tables, and a color-coded annotated mRNA sequence. Includes a copy-to-clipboard button for formatted sequence output. |
| `probes.tsv` | Tab-delimited probe table with oligo sequences ready for ordering. Columns: probe name, mRNA positions, CDS positions, length, target sequence, forward/reverse oligos, and mutations covered. |
| `{gene}_mutations.tsv` | Complete mutation data with computed mRNA positions, CDS positions, exon assignments, and probe coverage information. |

## Pipeline Details

### Data Sources

All data is fetched programmatically from NCBI E-utilities --- no manual downloads needed:

- **ClinVar mutations**: `esearch` + `esummary` for variants per gene (configurable molecular consequence filter)
- **Transcript structure**: Single GenBank `efetch` for the auto-discovered principal NM\_ transcript

An optional NCBI API key in `config.yaml` increases the rate limit from 3 to 10 requests/s. The fetcher includes automatic retry with exponential backoff on 429 (rate limit) errors, making multi-gene runs reliable even without an API key.

### Probe Design Algorithm

Greedy left-to-right sweep on the spliced mRNA:

1. Build forbidden zones around exon-exon junctions (configurable buffer)
2. Sort mutations by mRNA position
3. For each uncovered mutation, try all valid probe windows (target length, then min to max)
4. Score each window: coverage count x 100 + centering quality x 10 - overlap penalty x 50
5. Hard constraint: no overlap with forbidden zones
6. Soft constraint: overlap with existing probes is penalized, not rejected

### Oligo Format

| Oligo | Sequence |
|-------|----------|
| Forward (`{GENE}-{N}-FWD`) | `5'-AAAC-[reverse_complement(target)]-3'` |
| Reverse (`{GENE}-{N}-REV`) | `5'-AAAA-[target_sequence]-3'` |

### Coordinate System

- **mRNA position**: 1-based absolute position on the spliced mRNA (used for probe placement)
- **cDNA/CDS position** (`cdna_position`): the c.N value from HGVS notation, relative to the CDS start (c.1 = A of ATG start codon)
- **Conversion**: `mrna_position = cds_start + cdna_position - 1`

## Project Structure

```
pipeline/
  config.yaml              # Gene list and probe parameters
  Snakefile                # Master workflow definition
  rules/
    fetch_data.smk         # Data fetching rules (ClinVar, transcript)
    probe_design.smk       # Probe design rule
    report.smk             # Report generation and export rules
  scripts/
    fetch_clinvar.py       # NCBI ClinVar data fetcher
    fetch_transcript.py    # NCBI GenBank transcript fetcher
    design_probes.py       # Greedy probe placement algorithm
    generate_report.py     # Interactive HTML report generator
  envs/
    fetch.yaml             # Conda env: biopython, requests, pandas
    design.yaml            # Conda env: biopython, pandas
    report.yaml            # Conda env: jinja2, plotly, pandas
  data/                    # Intermediate files per gene (gitignored)
  results/                 # Final outputs per gene (gitignored)
  logs/                    # Per-rule log files
```

## Usage Examples

```bash
# Dry run (validate workflow without executing)
snakemake -n --cores 2

# Run full pipeline
snakemake --use-conda --cores 2

# Run up to a specific stage
snakemake data/POLR3A/clinvar_mutations.tsv --use-conda --cores 1
snakemake data/POLR3A/probes.tsv --use-conda --cores 1

# Generate a DAG visualization
snakemake --dag | dot -Tpng > dag.png
```

## Known Limitations

- ~9% of mutations cannot be covered (fall within junction buffer zones)
- Mutations at the very 5' end of the CDS (position < `scan_offset` from mRNA start) are skipped
- The probe overlap penalty is a heuristic; may need tuning for edge cases
- The HTML report uses Plotly CDN --- requires internet to load the JavaScript (data is embedded)
- Some ClinVar records may lack a parseable cDNA position --- these are kept in the mutations file but excluded from probe design
- Non-missense variants (e.g. nonsense) may lack a protein change --- these appear as empty strings in probe output

## License

This project is licensed under the [MIT License](LICENSE).

## Author

**Danny Bergeron** --- RNomics Platform
