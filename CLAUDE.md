# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Snakemake pipeline for designing antisense gene probes and oligos. Given a gene name in `config.yaml`, the pipeline automatically fetches ClinVar variants (configurable molecular consequence types) and transcript data from NCBI, designs probes (24-26 nt) that maximize mutation coverage while avoiding exon-exon junction buffer zones, and produces an interactive HTML report + orderable CSV.

Designed for any gene — POLR3A is the example, but changing `genes:` in config.yaml is all that's needed to target a different gene.

## Pipeline DAG

```
config.yaml
    │
    ▼
[fetch_clinvar] ──► data/{gene}/clinvar_mutations.tsv
    │                 data/{gene}/principal_transcript.txt
    ▼
[fetch_transcript] ──► data/{gene}/transcript.tsv  (includes cds_start, cds_end)
    │                    data/{gene}/junctions.tsv
    ▼
[design_probes] ──► data/{gene}/probes.tsv  (mrna_start/end + cds_start/end)
    │
    ├──► [generate_report] ──► results/{gene}/{gene}_report.html
    │                          results/{gene}/{gene}_mutations.tsv
    └──► [export_probes]     ──► results/{gene}/probes.tsv
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
    fetch_clinvar.py       # NCBI E-utilities: esearch + esummary for ClinVar variants (configurable consequence types)
    fetch_transcript.py    # NCBI efetch GenBank record → sequence + exon junctions + CDS extraction
    design_probes.py       # Greedy probe placement algorithm (converts cdna→mrna positions)
    generate_report.py     # Jinja2 + Plotly interactive HTML report + annotated sequence
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
      {gene}_mutations.tsv
      probes.tsv
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
snakemake results/POLR3A/POLR3A_mutations.tsv --use-conda --cores 1

# Generate DAG visualization
snakemake --dag | dot -Tpng > dag.png
```

## Config Schema (`config.yaml`)

```yaml
genes:
  - POLR3A          # List of gene symbols — pipeline runs independently per gene

ncbi_api_key: ""    # Optional: 10 req/s with key, 3 req/s without

molecular_consequences:               # ClinVar consequence filter — list, single value, or "all"
  - missense variant                   # Default if key omitted

probe_params:
  target_length: 25       # Ideal probe size (nt)
  min_length: 24           # Min acceptable
  max_length: 26           # Max acceptable
  junction_buffer: 5       # Excluded nt on each side of exon-exon junctions
  scan_offset: 5           # Start scanning this many nt from 5' end
```

The principal transcript is auto-discovered from ClinVar data (most common NM_ accession across all fetched variant records). No manual transcript ID needed.

## Data Sources (all NCBI E-utilities, fully programmatic)

1. **ClinVar mutations**: `esearch` + `esummary` for variants per gene
   - Molecular consequence filter is configurable via `molecular_consequences` in config.yaml
   - Accepts a single value (e.g. `["missense variant"]`), multiple values (e.g. `["missense variant", "nonsense"]`), or `"all"` (no consequence filter)
   - Multiple values are combined with OR in the esearch query
   - Valid ClinVar values: missense variant, intron variant, synonymous variant, non-coding transcript variant, frameshift variant, 5 prime UTR variant, 3 prime UTR variant, nonsense, splice donor variant, splice acceptor variant, inframe deletion, inframe insertion, initiator codon variant, inframe indel, stop lost, genic upstream transcript variant, genic downstream transcript variant, no sequence alteration
   - Returns: protein change, cDNA change, clinical significance, genomic position, ClinVar accession
   - cDNA position parsed from HGVS notation (`c.761G>A` → 761), stored as CDS-relative `cdna_position`
   - Batch esummary (200 IDs/call), rate limit 3 req/s without API key
   - Rate limiting: 0.5s delay between requests without API key (0.11s with key), with exponential backoff retry (up to 5 retries) on 429 errors — handles parallel multi-gene runs
   - Pagination supports >5000 results via `retstart`

2. **Transcript + exon + CDS structure**: Single GenBank `efetch` for NM_ accession
   - Returns mRNA sequence, exon features, and CDS feature with mRNA-relative positions
   - CDS start/end extracted from the GenBank CDS feature (e.g., CDS starts at mRNA position 109 for POLR3A)
   - Junctions computed from consecutive exon boundaries (e.g., exon 1 ends at 152, exon 2 starts at 153 → junction at position 152)

## Intermediate File Formats

### `data/{gene}/clinvar_mutations.tsv`
TSV with columns: `gene  transcript_id  protein_change  cdna_change  clinical_significance  chrom  pos_start  pos_end  clinvar_accession  cdna_position`
- `cdna_position` is the raw c.N value from HGVS cDNA notation (CDS-relative), NOT the absolute mRNA position
- Records without a parseable cDNA position have empty `cdna_position`

### `data/{gene}/principal_transcript.txt`
Single line: the auto-discovered NM_ accession (e.g., `NM_007055.4`)

### `data/{gene}/transcript.tsv`
Key-value TSV: `field\tvalue` with rows: `accession`, `sequence_length`, `organism`, `cds_start`, `cds_end`, `protein_start`, `sequence`
- `cds_start`/`cds_end`: 1-based mRNA positions of the CDS (coding sequence start/end)
- `protein_start`: first 30 amino acids of the protein translation (for verification)

### `data/{gene}/junctions.tsv`
TSV with columns: `junction_id  mrna_position  exon_before  exon_after`
- `mrna_position` is 1-based on the spliced mRNA

### `data/{gene}/probes.tsv`
TSV with columns: `probe_name  mrna_start  mrna_end  cds_start  cds_end  length  target_sequence  fwd_oligo  rev_oligo  mutations_covered  cdna_changes_covered  mutation_positions`
- `mrna_start`/`mrna_end`: absolute positions on the spliced mRNA (1-based)
- `cds_start`/`cds_end`: CDS-relative positions (can be negative for 5' UTR probes)
- `mutations_covered`: comma-separated protein changes (empty string if variant has no protein change, e.g. nonsense)
- `cdna_changes_covered`: comma-separated cDNA change strings (empty string if missing)
- `mutation_positions`: comma-separated mRNA positions
- `fwd_oligo` = `AAAC` + reverse_complement(target_sequence)
- `rev_oligo` = `AAAA` + target_sequence
- `probe_name` format: `{GENE}-{N}` (e.g., `POLR3A-1`)

## Coordinate System

- **mRNA position**: 1-based absolute position on the spliced mRNA transcript (used for probe placement, transcript map)
- **cDNA/CDS position** (`cdna_position`): the c.N value from HGVS notation, relative to the CDS start (c.1 = A of ATG start codon)
- **Conversion**: `mrna_position = cds_start + cdna_position - 1`
- `design_probes.py` reads `cds_start` from `transcript.tsv` and converts all positions before probe placement
- CDS positions can be negative for probes covering the 5' UTR

## Probe Design Algorithm

Greedy left-to-right sweep (all positions are mRNA-based):
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

Single self-contained file with Plotly CDN. Uses CSS classes throughout (no inline styles in visual HTML except dynamic probe background colors). Contains:
- **Summary**: gene, transcript, CDS range, mutation/probe counts, coverage %
- **Interactive transcript map**: 3-lane Plotly chart (exons, mutations, probes) with hover tooltips
- **ClinVar Mutations table**: 11 columns (ID, Protein Change, cDNA Change, Significance, mRNA Position, CDS Position, ClinVar Accession, Chrom, Genomic Pos, Exon, Probe), paginated, sortable, searchable, with horizontal scrollbar and CSV download
  - Mutation IDs use HGVS standard format: `NM_007055.4(POLR3A):c.251G>A` (transcript + gene in parentheses + cDNA change)
  - Probe column deduplicates — multiple mutations at the same position show the probe name once
- **Designed Probes table**: 11 columns (Probe Name, mRNA Start/End, CDS Start/End, Length, Target, FWD/REV Oligo, Mutations Covered, cDNA Changes), paginated (10/page default), sortable, searchable, with horizontal scrollbar and CSV download
  - cDNA Changes column has no max-width constraint (removed `seq-cell` class) so full content is always visible
- **Annotated mRNA Sequence**: exon blocks with labels, grey UTR backgrounds, colored mutation text with custom CSS tooltips (using `data-tip` attribute + `::after` pseudo-element for immediate display, `cursor: help`), probe regions shown as alternating colored backgrounds (`#dee6ef` for odd probes, `#ffffd7` for even probes, RGB-averaged blend `#eef2e3` where probes overlap), probe name labels centered at probe midpoints using invisible dash padding (`color: transparent; user-select: none`) in same 12px `'Courier New',monospace` font as nucleotides for exact character-level alignment, 100 nt per line in 10-nt groups
- **Copy Sequence button**: copies from a hidden div formatted for paste (60 nt/row, rich text with inline styles for paste compatibility, exon headers include probe ranges like "Exon 2 (positions 153-288, probes 2 to 5)", no position numbers or probe labels in sequence)

### Report CSS Classes

The visual HTML uses CSS classes (defined in the `<style>` block) instead of inline styles:
- **Nucleotides**: `.s` (base font), `.su` (UTR background), `.sm` (mutation bold), `.sm-p`/`.sm-lp`/`.sm-v`/`.sm-b`/`.sm-lb`/`.sm-o` (mutation color by significance), `.sg` (group separator)
- **Labels**: `.ll` (visible label text), `.lp` (invisible transparent dash padding, `user-select: none`)
- **Structure**: `.se` (exon container), `.sh` (exon header), `.sn` (nucleotide line), `.sl` (label line), `.sp` (empty spacer)
- **Template**: `.summary-info`, `.table-wrap`, `.copy-hidden`, `.ld-*` (legend dots), `.slb-*` (seq legend boxes), `.slm` (seq legend mutation)
- Probe backgrounds stay inline (`style="background:..."`) since they're dynamic hex values
- The copy HTML retains full inline styles for rich text paste compatibility

## Output Files

- `results/{gene}/{gene}_report.html`: Interactive HTML report
- `results/{gene}/probes.tsv`: Probe data (same columns as report probes table)
- `results/{gene}/{gene}_mutations.tsv`: Detailed mutation data with all computed fields (same as CSV download from report)

## Known Limitations / Improvement Areas

- ~9% of mutations cannot be covered (fall within junction buffer zones)
- Mutations at the very 5' end of the CDS (position < scan_offset from mRNA start) are skipped
- The probe overlap penalty (−50) is a heuristic; may need tuning for edge cases
- Report uses Plotly CDN — requires internet to load JS (data is embedded)
- Some ClinVar records may lack a parseable cDNA position — these are kept in the mutations file but excluded from probe design
- Non-missense variants (e.g. nonsense) may lack a protein change — these are stored as empty strings in probe output
- The Copy Sequence button uses Clipboard API which may not work in all browsers (falls back to execCommand)

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
