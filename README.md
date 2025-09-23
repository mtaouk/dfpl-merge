# dfpl-merge

`dfpl-merge` is a command-line tool for merging, consolidating, and
resolving DefenseFinder and PADLOC results with genomic coordinates from
Bakta annotations which helps for visualization and other exploratory
work. It produces both a full merged table and a consolidated summary
table, simplifying downstream analyses of defense systems in prokaryotic
genomes.

------------------------------------------------------------------------

## Features

-   Merges **DefenseFinder** and **PADLOC** outputs on `locus_tag`.
-   Fills in missing genomic coordinates using **Bakta CDS annotations**
    (if provided). Note: input must be in the **Bakta** tsv format, but
    **Prokka** annotations can be input if first transformed to the
    correct format (see below).
-   Consolidates gene and system names into a single, harmonised set.
-   Produces both a full merged table and a summary table of source and
    consolidated information.

------------------------------------------------------------------------

## Installation

This tool is written in Python 3 and relies on the following packages:

-   `pandas` (for dataframe manipulation)
-   Standard library modules: `argparse`, `csv`, `pathlib`

### Steps

1.  Clone this repository:

``` bash
git clone https://github.com/mtaouk/dfpl-merge.git
cd dfpl-merge
```

2.  Install the required package

``` bash
conda install pandas
```

or

``` bash
pip install pandas
```

3.  Run the script

------------------------------------------------------------------------

## Usage

``` bash
options:
  -d DEFENSEFINDER_TSV  DefenseFinder genes table (required)
  -p PADLOC_TSV         PADLOC results table (required)
  -o OUTDIR             Output directory (required)
  -b BAKTA_TSV          Bakta annotations (optional)
  -v, --version         show program's version number and exit
  -h, --help            Show this help and exit
```

### Example

``` bash
python dfpl-merge.py \
  -d path/to/defensefinder.tsv \
  -p path/to/padloc.tsv \
  -b path/to/bakta.tsv \
  -o path/to/output_dir \
```

------------------------------------------------------------------------

## Inputs

The tool accepts three inputs:

1.  **DefenseFinder TSV** (`-d`)
    -   defense_finder_genes.tsv output tsv from DefenseFinder.
    -   Must include columns: `hit_id`, `gene_name`, `sys_id`,
        `hit_i_eval`, `hit_profile_cov`, `hit_seq_cov`, `model_fqn`,
        `hit_status`, `sys_wholeness`, `hit_score`.
    -   Recommended usage of Defense-finder:
        `defense-finder run annotated_genome.faa -o output_dir`
2.  **PADLOC CSV** (`-p`)
    -   Standard output csv from PADLOC.
    -   Must include columns: `system`, `target.name`, `target.name`,
        `protein.name`, `full.seq.E.value`, `domain.iE.value`,
        `target.coverage`, `hmm.coverage`, `start`, `end`, `strand`.
    -   Recommended usage of Defense-finder:
        `padloc --faa annotated_genome.faa --gff annotated_genome.gff3`
3.  **Bakta annotation TSV** (`-b`) (Optional)
    -   Genomic annotations with coding sequences.
    -   Must contain a commented header line beginning with
        `#Sequence ...` and columns: `Locus Tag`, `Start`, `End`,
        `Strand`.
    -   If this is not provided, it just means that any genes identified
        by DefenseFinder and not PADLOC will not have any coordinates
        filled in for them. Prokka annotations can be input if they are
        restructured into this format first.

------------------------------------------------------------------------

## Important

-   This tool **only works if you have used the same assembly as input
    for both** `DefenseFinder` **and** `PADLOC`.
-   `Bakta` annotated assemblies are preferred.

------------------------------------------------------------------------

## Logic / Workflow

1.  **Load and clean inputs:**
    -   DefenseFinder: rename columns, extract gene/system names,
        simplify model identifiers.
    -   PADLOC: rename columns, normalise strand, cast coordinates to
        integers.
    -   Bakta: parse CDS records, retain only essential coordinates. (If
        provided).
2.  **Merge tables:**
    -   Outer merge on `locus_tag` to include all entries from both
        DefenseFinder and PADLOC.
    -   Coordinates from PADLOC are used if present; otherwise, Bakta
        coordinates backfill missing values. If not Bakta table is
        provided, genes that are identified with DefenseFinder only will
        have missing values for `Start`, `End`, and `Strand`.
3.  **Generate outputs:**
    -   **Full merged table:** key columns from DefenseFinder and
        PADLOC, plus unified `start`, `end`, `strand`.
    -   **Consolidated summary table:**
        -   `locus_tag`
        -   `source_type` → `"DefenseFinder only"`, `"PADLOC only"`, or
            `"both"`
        -   `consolidated_gene` → DefenseFinder gene name if present,
            otherwise PADLOC gene name
        -   `consolidated_system` → DefenseFinder system name if
            present, otherwise PADLOC system

------------------------------------------------------------------------

## Outputs

### 1. Full Merged Table (`defensefinder_padloc_merged.tsv`)

| Column                          | Description                        |
|---------------------------------|------------------------------------|
| `locus_tag`                     | Unique gene identifier             |
| `defensefinder_model`           | DefenseFinder model name           |
| `defensefinder_gene_name`       | DefenseFinder gene name            |
| `defensefinder_system`          | DefenseFinder system identifier    |
| `defensefinder_hit_i_eval`      | DefenseFinder e-value for hit      |
| `defensefinder_hit_profile_cov` | Profile coverage                   |
| `defensefinder_hit_seq_cov`     | Sequence coverage                  |
| `defensefinder_hit_status`      | Hit status                         |
| `defensefinder_sys_wholeness`   | System completeness                |
| `defensefinder_hit_score`       | DefenseFinder hit score            |
| `padloc_system`                 | PADLOC system identifier           |
| `padloc_gene_name`              | PADLOC gene name                   |
| `padloc_evalue`                 | PADLOC full sequence e-value       |
| `padloc_domain_ievalue`         | PADLOC domain e-value              |
| `padloc_target_cov`             | PADLOC target coverage             |
| `padloc_hmm_cov`                | PADLOC HMM coverage                |
| `start`                         | Start coordinate (PADLOC or Bakta) |
| `end`                           | End coordinate (PADLOC or Bakta)   |
| `strand`                        | Strand (+ / -)                     |

### 2. Consolidated Summary Table (`defensefinder_padloc_consolidated.tsv`)

| Column                | Description                                             |
|----------------------------|--------------------------------------------|
| `locus_tag`           | Unique gene identifier                                  |
| `source_type`         | `"DefenseFinder only"`, `"PADLOC only"`, or `"both"`    |
| `consolidated_gene`   | Unified gene name (prefer DefenseFinder if available)   |
| `consolidated_system` | Unified system name (prefer DefenseFinder if available) |

------------------------------------------------------------------------
