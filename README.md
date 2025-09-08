# df-padloc-merge

`df-padloc-merge` is a command-line tool for merging, consolidating, and resolving DefenseFinder and PADLOC results with genomic coordinates from Bakta annotations. It produces both a full merged table and a consolidated summary table, simplifying downstream analyses of defence systems in prokaryotic genomes.

---

## Features

- Merges **DefenseFinder** and **PADLOC** outputs on `locus_tag`.
- Fills in missing genomic coordinates using **Bakta CDS annotations**.
- Consolidates gene and system names into a single, harmonised set.
- Produces both a full merged table and a summary table of source and consolidated information.
- Handles strand normalisation and coordinate type casting automatically.

---

## Inputs

The tool requires three inputs:

1. **DefenseFinder TSV** (`-d` / `--defensefinder`)  
   - defense_finder_genes.tsv output tsv from DefenseFinder.  
   - Must include columns: `replicon`, `hit_id`, `gene_name`, `sys_id`, `hit_i_eval`, `hit_profile_cov`, `hit_seq_cov`, `model_fqn`, `hit_status`, `sys_wholeness`, `hit_score`.  

2. **PADLOC CSV** (`-p` / `--padloc`)  
   - Standard output csv from PADLOC.  
   - Must include columns: `system`, `target.name`, `start`, `end`, `strand`.  

3. **Bakta annotation TSV** (`-b` / `--bakta`)  
   - Genomic annotations with coding sequences.  
   - Must contain a commented header line beginning with `#Sequence ...` and columns: `Locus Tag`, `Start`, `Stop`, `Strand`.  

---

## Important

This tool **only works if you have used Bakta to annotate your assemblies** and then run DefenseFinder and PADLOC on that same assembly.  
It has to be the **same genome** or it won't work. If the genome wasn't annotated with Bakta, coordinates cannot be filled, and the merge will fail.

---

## Logic / Workflow

1. **Load and clean inputs:**  
   - DefenseFinder: rename columns, extract gene/system names, simplify model identifiers.  
   - PADLOC: rename columns, normalise strand, cast coordinates to integers.  
   - Bakta: parse CDS records, retain only essential coordinates.  

2. **Merge tables:**  
   - Outer merge on `locus_tag` to include all entries from both DefenseFinder and PADLOC.  
   - Coordinates from PADLOC are used if present; otherwise, Bakta coordinates backfill missing values.  

3. **Generate outputs:**  
   - **Full merged table:** key columns from DefenseFinder and PADLOC, plus unified `start`, `end`, `strand`.  
   - **Consolidated summary table:**  
     - `locus_tag`  
     - `source_type` → `"DF only"`, `"PADLOC only"`, or `"both"`  
     - `consolidated_gene` → DefenseFinder gene name if present, otherwise PADLOC gene name  
     - `consolidated_system` → DefenseFinder system name if present, otherwise PADLOC system

---

## Outputs

### 1. Full Merged Table (`df_padloc_merged.tsv`)

| Column | Description |
|--------|-------------|
| `locus_tag` | Unique gene identifier |
| `df_model` | DefenseFinder model name |
| `df_gene_name` | DefenseFinder gene name |
| `df_system` | DefenseFinder system identifier |
| `df_hit_i_eval` | DefenseFinder e-value for hit |
| `df_hit_profile_cov` | Profile coverage |
| `df_hit_seq_cov` | Sequence coverage |
| `df_hit_status` | Hit status |
| `df_sys_wholeness` | System completeness |
| `df_hit_score` | DefenseFinder hit score |
| `sample_name` | Replicon/sample identifier |
| `start` | Start coordinate (PADLOC or Bakta) |
| `end` | End coordinate (PADLOC or Bakta) |
| `strand` | Strand (+ / -) |
| `pl_system` | PADLOC system identifier |
| `pl_gene_name` | PADLOC gene name |
| `pl_evalue` | PADLOC full sequence e-value |
| `pl_domain_ievalue` | PADLOC domain e-value |
| `pl_target_cov` | PADLOC target coverage |
| `pl_hmm_cov` | PADLOC HMM coverage |

### 2. Consolidated Summary Table (`df_padloc_consolidated.tsv`)

| Column | Description |
|--------|-------------|
| `locus_tag` | Unique gene identifier |
| `source_type` | `"DF only"`, `"PADLOC only"`, or `"both"` |
| `consolidated_gene` | Unified gene name (prefer DF if available) |
| `consolidated_system` | Unified system name (prefer DF if available) |

---

## Usage

```bash
python df-padloc-merge.py \
  -d path/to/defensefinder.tsv \
  -p path/to/padloc.tsv \
  -b path/to/bakta.tsv \
  -o path/to/output_dir \
  -v

