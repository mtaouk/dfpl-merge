#!/usr/bin/env python3
 
"""
Copyright 2025 Mona Taouk
https://github.com/mtaouk/dfpl-merge.py

This is the main script for dfpl-merge. 

dfpl-merge is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. dfpl-merge is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with dfpl-merge. If
not, see <http://www.gnu.org/licenses/>.
"""

# ---- imports ---- 
import argparse # parse command-line flags
import csv # for reading the input files
from pathlib import Path # nicer file paths than raw strings 
import sys # for writing errors
import pandas as pd # aparently need this for the merging of dataframes

# ---- main ----
def main():
    args = build_parser().parse_args()

    # load + tidy
    defensefinder_rows = read_tidy_defensefinder(args.d)
    padloc_rows = read_tidy_padloc(args.p)
    bakta_rows = read_tidy_bakta(args.b) if args.b else None

    merged = merge_all(defensefinder_rows, padloc_rows, bakta_records=bakta_rows)

    # after merging, format values
    for r in merged:
        # format e-values nicely
        for col in ("padloc_domain_ievalue", "padloc_evalue"):
            val = r.get(col)
            if val not in (None, "NA", ""):
                try:
                    r[col] = f"{float(val):.2e}"
                except ValueError:
                    pass

        # force start/end to full integer strings
        for col in ("start", "end"):
            val = r.get(col)
            if val not in (None, "NA", ""):
                try:
                    r[col] = str(int(float(val)))
                except ValueError:
                    pass

    outpath = args.o / "defensefinder_padloc_merged.tsv"
    
    # define field order with coordinates last
    if merged:
        cols = [c for c in merged[0].keys() if c not in ("start", "end", "strand")]
        field_order = cols + ["start", "end", "strand"]
    else:
        field_order = ["locus_tag", "defensefinder_model", "defensefinder_gene_name", "defensefinder_system", "defensefinder_hit_i_eval", 
                        "defensefinder_hit_profile_cov", "defensefinder_hit_seq_cov", "defensefinder_hit_status", "defensefinder_sys_wholeness",
                        "defensefinder_hit_score", "sample_name", "padloc_system", "padloc_gene_name", "padloc_evalue",
                        "padloc_domain_ievalue", "padloc_target_cov", "padloc_hmm_cov", "start", "end", "strand"]
    write_tsv(merged, outpath, field_order=field_order)

    summary = make_summary_table(merged)
    summary_path = args.o / "defensefinder_padloc_consolidated.tsv"
    write_tsv(summary, summary_path, field_order=["locus_tag", "source_type", "consolidated_gene", "consolidated_system"])
    

# ----argument parser ---- 
def build_parser():
    p = argparse.ArgumentParser(
        prog="df-padloc-merge",
        description="Merge, consolidate and resolve DefenseFinder and Padloc outputs",
        epilog="Example: dfpl-merge.py -d defensefinder.tsv -p PADLOC.tsv -b bakta.tsv -o out",
        add_help=False
    )

    # required arguments
    p.add_argument("-d", metavar="DEFENSEFINDER_TSV", type=Path, required=True, help="DefenseFinder genes table (required)")
    p.add_argument("-p", metavar="PADLOC_TSV", type=Path, required=True, help="PADLOC results table (required)")
    p.add_argument("-o", metavar="OUTDIR", type=Path, required=True, help="Output directory (required)")

    # optional argument
    p.add_argument("-b", metavar="BAKTA_TSV", type=Path, required=False, help="Bakta annotations (optional)")
    # version flag that works standalone
    p.add_argument("-v", "--version", action="version", version="v0.1beta")
    p.add_argument("-h", "--help", action="help", help="Show this help and exit")

    return p


# ---- generic table reader ----
def read_table(path, required=None, keep=None):
    """Read CSV or TSV into list[dict]. Validate minimal columns. Optionally keep a subset."""
    path = Path(path)
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open("r", newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter=delim)
        if required:
            missing = [c for c in required if c not in rdr.fieldnames]
            if missing:
                raise ValueError(f"{path} missing required columns: {missing}\nFound: {rdr.fieldnames}")
        rows = list(rdr)
    if keep:
        rows = [{k: r.get(k) for k in keep if k in r} for r in rows]
    return rows


# ---- clean Defensefinder columns ----
def read_tidy_defensefinder(path, keep_only_mapped=True):
    keep = ["hit_id", "gene_name", "sys_id", "hit_i_eval", "hit_profile_cov", "hit_seq_cov", "model_fqn", "hit_status", "sys_wholeness", "hit_score"]
    rows = read_table(path, keep, keep)
    if not rows:
        sys.exit("Error: DEFENSEFINDER_TSV file is empty")

    ren = {
        "hit_id": "locus_tag",
        "model_fqn": "defensefinder_model",
        "gene_name": "defensefinder_gene_name",
        "sys_id": "defensefinder_system",
        "hit_i_eval": "defensefinder_hit_i_eval",
        "hit_profile_cov": "defensefinder_hit_profile_cov",
        "hit_seq_cov": "defensefinder_hit_seq_cov",
        "hit_status": "defensefinder_hit_status",
        "sys_wholeness": "defensefinder_sys_wholeness",
        "hit_score": "definder_hit_score"
    }

    out = []
    for rec in rows:
        r = rec.copy()

        # gene_name
        gn = r.get("gene_name", "")
        if "__" in gn:
            r["gene_name"] = gn.split("__", 1)[1]

        # derive from model_fqn
        model = r.get("model_fqn") or ""
        if model:
            # sys_id = last path token
            last_tok = model.strip("/").split("/")[-1]
            if last_tok:
                r["sys_id"] = last_tok
            # model_fqn = token after 'defense-finder-models/'
            parts = model.strip("/").split("/")
            try:
                i = parts.index("defense-finder-models")
                if i + 1 < len(parts):
                    r["model_fqn"] = parts[i + 1]
            except ValueError:
                pass  # leave as-is if anchor missing

        # rename keys
        if keep_only_mapped:
            r = {ren[k]: r[k] for k in ren if k in r}
        else:
            for old, new in ren.items():
                if old in r:
                    r[new] = r.pop(old)

        out.append(r)

    return out

# ---- clean Padloc columns ----
def read_tidy_padloc(path, keep_only_mapped=True):
    required = ["system","target.name","start","end","strand"]
    keep = required + ["protein.name", "full.seq.E.value","domain.iE.value","target.coverage","hmm.coverage"]
    rows = read_table(path, required, keep)
    if not rows:
        sys.exit("Error: PADLOC_TSV file is empty")

    ren = {
        "target.name": "locus_tag",
        "system": "padloc_system",
        "protein.name": "padloc_gene_name",
        "full.seq.E.value": "padloc_evalue",
        "domain.iE.value": "padloc_domain_ievalue",
        "target.coverage": "padloc_target_cov",
        "hmm.coverage": "padloc_hmm_cov",
        "start": "start",
        "end": "end",
        "strand": "strand"
    }

    out = []
    for rec in rows:
        r = rec.copy()

        # rename keys
        if keep_only_mapped:
            r = {new: r.get(old, "") for old, new in ren.items() if old in r}
        else:
            for old, new in ren.items():
                if old in r and new != old:
                    r[new] = r.pop(old)

        # cast coordinates
        for k in ("start", "end"):
            if k in r:
                try:
                    # sometimes PADLOC writes floats-as-strings; be chill
                    r[k] = int(float(r[k]))
                except Exception:
                    pass

        # normalise strand
        if "strand" in r:
            s = str(r["strand"]).strip()
            if s in ("1", "+", "plus"):
                r["strand"] = "+"
            elif s in ("-1", "-", "minus"):
                r["strand"] = "-"
            # else leave whatever PADLOC gave

        out.append(r)
    return out

# ---- Bakta TSV reader and cleaner ----
def read_tidy_bakta(path):
    """
    Read a Bakta 'Annotated with Bakta' TSV (with a commented header),
    keep only selected columns, rename them, and return CDS rows as a list of dicts.
    """
    # only keep/rename these columns
    ren = {
        "Locus Tag": "locus_tag",
        "Start": "start",
        "Stop": "end",
        "Strand": "strand"
    }
    path = Path(path)
    with path.open("r", newline="", encoding="utf-8") as fh:
        header = None
        # find the commented header line (e.g. "#Sequence Id\tType\tStart...")
        for line in fh:
            if not line.strip():
                continue
            text = line.lstrip()
            if text.startswith("#Sequence"):
                header = [h.strip() for h in text.lstrip("#").strip().split("\t")]
                break
        if header is None:
            raise ValueError(f"{path} has no '#Sequence ...' header line")
        # make sure all required cols exist
        missing = [k for k in ren if k not in header]
        if missing:
            raise ValueError(f"{path} missing required columns: {missing}")
        # now read only data lines (skip other comment lines)
        data_iter = (ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#"))
        reader = csv.DictReader(data_iter, delimiter="\t", fieldnames=header)
        out = []
        for row in reader:
            if row.get("Type") != "cds":
                continue  # keep only coding sequences
            # keep only columns in ren and rename them
            rec = {ren[src]: row[src] for src in ren}
            out.append(rec)
    if not out:
        sys.exit("Error: BAKTA_TSV file is empty")
    return out

# ---- merge DF + Padloc (+ Bakta coords) ----
def merge_all(defensefinder_rows, padloc_rows, bakta_records=None):
    defensefinder = pd.DataFrame(defensefinder_rows)
    padloc = pd.DataFrame(padloc_rows)

    # outer merge on locus_tag 
    merged = pd.merge(defensefinder, padloc, on="locus_tag", how="outer")

    # backfill with Bakta coords if provided
    if bakta_records is not None:
        bakta = pd.DataFrame(bakta_records)
        merged = pd.merge(merged, bakta, on="locus_tag", how="left", suffixes=("", "_bakta"))

        for col in ["start", "end", "strand"]:
            bakta_col = f"{col}_bakta"
            if bakta_col in merged.columns:
                merged[col] = merged[col].where(merged[col].notna(), merged[bakta_col])

        # drop helper bakta columns
        merged = merged.drop(columns=[c for c in merged.columns if c.endswith("_bakta")])

    # replace NaN with "NA" for output
    merged = merged.fillna("NA")

    return merged.to_dict(orient="records")

# ---- summary table ----
def make_summary_table(merged_records):
    df = pd.DataFrame(merged_records)

    def source_type(row):
        has_defensefinder = row.get("defensefinder_gene_name") not in (None, "NA", "")
        has_padloc = row.get("padloc_gene_name") not in (None, "NA", "")
        if has_defensefinder and has_padloc:
            return "Both"
        elif has_defensefinder:
            return "DefenseFinder only"
        elif has_padloc:
            return "Padloc only"
        else:
            return "NA"

    summary = pd.DataFrame()
    summary["locus_tag"] = df["locus_tag"]
    summary["source_type"] = df.apply(source_type, axis=1)
    summary["consolidated_gene"] = df.apply(
        lambda r: r["defensefinder_gene_name"] if r["defensefinder_gene_name"] not in (None, "NA", "") else r["padloc_gene_name"],
        axis=1,
    )
    summary["consolidated_system"] = df.apply(
        lambda r: r["defensefinder_system"] if r["defensefinder_system"] not in (None, "NA", "") else r["padloc_system"],
        axis=1,
    )

    return summary.to_dict(orient="records")


# ---- write TSVs ----
def write_tsv(rows, path, field_order):
    """Always write a header using field_order, then rows (if any)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=list(field_order),
            delimiter="\t",
            extrasaction="ignore"  # ignore any unexpected keys in rows
        )
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in field_order})

if __name__ == "__main__":
    main()