# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Snakemake pipeline for designing antisense gene probes and oligos. Given a gene name in `config.yaml`, the pipeline automatically fetches ClinVar missense mutations and transcript data from NCBI, designs probes (24-26 nt) that maximize mutation coverage while avoiding exon-exon junction buffer zones, and produces an interactive HTML report + orderable CSV.

Designed for any gene — POLR3A is the example, but changing `genes:` in config.yaml is all that's needed to target a different gene.

## Pipeline DAG

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
    └──► [export_probes]     ──► results/{gene}/probes.csv
```

## File Structure

```
pipeline/
  config.yaml              # Gene list, probe params, optional NCBI API key
  Snakefile                # Master workflow: configfile + includes + rule all
  rules/
    fetch_data.smk         # fetch_clinvar, fetch_transcript
    probe_design.smk       # design_probes
    report.smk             # generate_report, export_probes
  scripts/
    fetch_clinvar.py       # NCBI E-utilities: esearch + esummary for ClinVar missense variants
    fetch_transcript.py    # NCBI efetch GenBank record → sequence + exon junction extraction
    design_probes.py       # Greedy probe placement algorithm
    generate_report.py     # Jinja2 + Plotly interactive HTML report
  envs/
    fetch.yaml             # biopython, requests, pandas
    design.yaml            # biopython, pandas
    report.yaml            # jinja2, plotly, pandas
  data/                    # Intermediate files per gene (gitignored)
    {gene}/
      clinvar_mutations.tsv
      principal_transcript.txt
      transcript.tsv
      junctions.tsv
      probes.tsv
  results/                 # Final outputs per gene (gitignored)
    {gene}/
      {gene}_report.html
      probes.csv
  logs/                    # Per-rule log files
```

## Environment Setup

Before running any snakemake command, source the shell config and activate the snakemake conda environment:
```bash
source /home/danx/.zshrc
conda activate snakemake7.24
```

## Commands

All snakemake commands require `--cores <N>` to work properly.

```bash
# Dry run (validate workflow without executing)
snakemake -n --cores 2

# Run full pipeline
snakemake --use-conda --cores 2

# Run a specific output (rules use wildcards, target files instead of rule names)
snakemake data/POLR3A/clinvar_mutations.tsv --use-conda --cores 1
snakemake data/POLR3A/probes.tsv --use-conda --cores 1
snakemake results/POLR3A/POLR3A_report.html --use-conda --cores 1

# Generate DAG visualization
snakemake --dag | dot -Tpng > dag.png
```

## Config Schema (`config.yaml`)

```yaml
genes:
  - POLR3A          # List of gene symbols — pipeline runs independently per gene

ncbi_api_key: ""    # Optional: 10 req/s with key, 3 req/s without

probe_params:
  target_length: 25       # Ideal probe size (nt)
  min_length: 24           # Min acceptable
  max_length: 26           # Max acceptable
  junction_buffer: 5       # Excluded nt on each side of exon-exon junctions
  scan_offset: 5           # Start scanning this many nt from 5' end
```

The principal transcript is auto-discovered from ClinVar data (most common NM_ accession across all missense records). No manual transcript ID needed.

## Data Sources (all NCBI E-utilities, fully programmatic)

1. **ClinVar mutations**: `esearch` + `esummary` for missense variants
   - Search term: `"{gene}"[Gene] AND single_gene[prop] AND "missense variant"[molecular_consequence]`
   - Returns: protein change, cDNA change, clinical significance, genomic position, ClinVar accession
   - mRNA position parsed from cDNA notation (`c.761G>A` → 761)
   - Batch esummary (200 IDs/call), rate limit 3 req/s without API key

2. **Transcript + exon structure**: Single GenBank `efetch` for NM_ accession
   - Returns mRNA sequence AND exon features with mRNA-relative positions
   - Junctions computed from consecutive exon boundaries (e.g., exon 1 ends at 152, exon 2 starts at 153 → junction at position 152)

## Intermediate File Formats

### `data/{gene}/clinvar_mutations.tsv`
TSV with columns: `gene  transcript_id  protein_change  cdna_change  clinical_significance  chrom  pos_start  pos_end  clinvar_accession  mrna_position`

### `data/{gene}/principal_transcript.txt`
Single line: the auto-discovered NM_ accession (e.g., `NM_007055.4`)

### `data/{gene}/transcript.tsv`
Key-value TSV: `field\tvalue` with rows: `accession`, `sequence_length`, `organism`, `sequence`

### `data/{gene}/junctions.tsv`
TSV with columns: `junction_id  mrna_position  exon_before  exon_after`
- `mrna_position` is 1-based on the spliced mRNA

### `data/{gene}/probes.tsv`
TSV with columns: `probe_name  start  end  length  target_sequence  fwd_oligo  rev_oligo  mutations_covered  mutation_positions`
- `fwd_oligo` = `AAAC` + reverse_complement(target_sequence)
- `rev_oligo` = `AAAA` + target_sequence
- `probe_name` format: `{GENE}-{N}` (e.g., `POLR3A-1`)
- Mutations in semicolon-separated format

## Probe Design Algorithm

Greedy left-to-right sweep:
1. Build forbidden zones: junction_buffer nt on each side of every exon-exon junction
2. Sort mutations by mRNA position
3. For each uncovered mutation:
   - Try all valid windows (target_length, then min to max) containing the mutation
   - Score: coverage count × 100 + centering quality × 10 − overlap penalty × 50
   - Centering measured against midpoint of covered mutation range (first + last) / 2, not just the first mutation
   - Hard constraint: no overlap with forbidden zones
   - Soft constraint: no overlap with existing probes (penalized, not rejected)
4. Unplaceable mutations (in forbidden zones or at edges) logged as warnings

## Probe Output Naming

- Forward oligo: `{GENE}-{N}-FWD` → `5'-AAAC-[reverse_complement(target)]-3'`
- Reverse oligo: `{GENE}-{N}-REV` → `5'-AAAA-[target_sequence]-3'`

## HTML Report

Single self-contained file with Plotly CDN. Contains:
- **Summary**: gene, transcript, mutation/probe counts, coverage %
- **Interactive transcript map**: 3-lane Plotly chart (exons, mutations, probes) with hover tooltips
- **Mutation table**: paginated (10/page), sortable, searchable, ClinVar details
- **Probe table**: with CSV download button

## Known Limitations / Improvement Areas

- ~9% of mutations cannot be covered (fall within junction buffer zones)
- The CSV download in the report generates FWD/REV rows with `{probe_name}-FWD` and `{probe_name}-REV` naming
- Mutations at the very 5' end (position < scan_offset) are skipped
- The probe overlap penalty (−50) is a heuristic; may need tuning for edge cases
- Report uses Plotly CDN — requires internet to load JS (data is embedded)

## Reference Documents

Project specifications are in `../infos/`:
- `POLR3A-ClinVar-Missense MUTATIONS.xlsx` — mutation annotations
- `POLR3A_Target Sequences.docx` — target sequences
- `TARGET SEQUENCES-OLIGOS to ORDER.xlsx` — oligo specifications
- `WORKFLOW for the design of gRNAs.docx` — methodology reference
- `POLR3A_annotated Mutations.docx` — annotated mutation catalog
- `POLR3A-ccDS Sequence.docx` — coding sequence reference

## Snakemake Conventions

- Rules in `rules/*.smk` are grouped by concern, NOT one file per rule
- Conda env paths in rule files must be `../envs/*.yaml` (resolved relative to `rules/`)
- Script paths in rule files must be `../scripts/*.py` (same reason)
- `log:` directive must come before `conda:` and `script:` in rule definitions
- Rules with wildcards cannot be called by name — target their output files instead
- The `run:` directive cannot use `conda:` — use `shell:` instead

## Conda Channel Priority

Conda may warn about strict channel priority. This is cosmetic and does not affect pipeline runs. If desired, fix with:
```bash
conda config --set channel_priority strict
```
