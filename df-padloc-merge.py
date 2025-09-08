#!/usr/bin/env python3


# ---- imports ---- 
import argparse # parse command-line flags
import csv # for reading the input files
from pathlib import Path # nicer file paths than raw strings 
import pandas as pd # aparently need this for the merging of dataframes


# ---- pretty help with single metavar ---- 
class OneMetavarHelp(argparse.HelpFormatter):
    def __init__(self, prog):
        super().__init__(prog, max_help_position=36, width=100)
    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        opts = ", ".join(action.option_strings)
        if action.nargs == 0:
            return opts
        mv = self._metavar_formatter(action, action.dest)(1)[0]
        return f"{opts} {mv}"


# ----argument parser ---- 
def build_parser():
    p = argparse.ArgumentParser(
        prog="df-padloc-merge", 
        description="Merge, consolidate and resolve DefenseFinder and Padloc outputs",
        epilog="Example: df-padloc-merge -d DF.tsv -p PADLOC.tsv -b bakta.tsv -o out",
        add_help=False,
        formatter_class=OneMetavarHelp)
    
    req = p.add_argument_group("Required arguments")
    req.add_argument("-d", "--defensefinder", metavar="DF_TSV", type=Path, required=True, help="DefenseFinder genes table")
    req.add_argument("-p", "--padloc", metavar="PADLOC_TSV", type=Path, required=True, help="PADLOC results table")
    req.add_argument("-b", "--bakta", metavar="BAKTA_GFF_OR_TSV", type=Path, required=True, help="Bakta annotations")
    req.add_argument("-o", "--outdir", metavar="OUTDIR", type=Path, required=True, help="Output directory")
   
    opt = p.add_argument_group("Optional arguments")
    opt.add_argument("-h", "--help", action="help", help="Show this help and exit")
    opt.add_argument("-v", "--verbose", action="store_true", help="Verbose logs")
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


# ---- specific loaders (thin wrappers) ----
def load_df_genes(path):
    keep = ["replicon", "hit_id", "gene_name", "sys_id", "hit_i_eval", "hit_profile_cov", "hit_seq_cov", "model_fqn", "hit_status", "sys_wholeness", "hit_score"]
    return read_table(path, required=keep, keep=keep)

def load_padloc(path):
    required = ["system","target.name","start","end","strand"]
    keep = required + ["protein.name", "full.seq.E.value","domain.iE.value","target.coverage","hmm.coverage"]
    return read_table(path, required=required, keep=keep)


# ---- clean Defensefinder columns ----
def tidy_df_genes(rows, keep_only_mapped=False):
    """
    Clean DefenceFinder gene rows in memory.
    - gene_name: keep suffix after "__"
    - sys_id: replace with last token of model_fqn path
    - model_fqn: reduce to the token after 'defense-finder-models/'
    - rename keys per `ren`, keeping others by default
    """
    ren = {
        "hit_id": "locus_tag",
        "model_fqn": "df_model",
        "gene_name": "df_gene_name",
        "sys_id": "df_system",
        "hit_i_eval": "df_hit_i_eval",
        "hit_profile_cov": "df_hit_profile_cov",
        "hit_seq_cov": "df_hit_seq_cov",
        "hit_status": "df_hit_status",
        "sys_wholeness": "df_sys_wholeness",
        "hit_score": "df_hit_score"
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
def tidy_padloc(rows, keep_only_mapped=False):
    """
    Clean PADLOC rows in memory:
      - rename PADLOC headers to unified names
      - cast start/end to int where possible
      - normalise strand to '+' / '-'
    """
    ren = {
        "target.name": "locus_tag",
        "system": "pl_system",
        "protein.name": "pl_gene_name",
        "full.seq.E.value": "pl_evalue",
        "domain.iE.value": "pl_domain_ievalue",
        "target.coverage": "pl_target_cov",
        "hmm.coverage": "pl_hmm_cov",
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
def load_bakta_clean(path):
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
            # cast coords
            try:
                rec["start"] = int(rec["start"])
                rec["end"] = int(rec["end"])
            except Exception:
                pass
            out.append(rec)

    return out

# ---- merge DF + Padloc (+ Bakta coords) ----
def merge_all(df_rows, pl_rows, bakta_records=None):
    df = pd.DataFrame(df_rows)
    pl = pd.DataFrame(pl_rows)

    # outer merge on locus_tag
    merged = pd.merge(df, pl, on="locus_tag", how="outer", suffixes=("_df", "_pl"))

    # backfill with Bakta coords if provided
    if bakta_records is not None:
        bakta = pd.DataFrame(bakta_records)
        merged = pd.merge(merged, bakta, on="locus_tag", how="left", suffixes=("", "_bakta"))

        for col in ["start", "end", "strand"]:
            bakta_col = f"{col}_bakta"
            if bakta_col in merged.columns:
                merged[col] = merged.get(col).where(merged.get(col).notna(), merged[bakta_col])

        # drop helper bakta columns
        merged = merged.drop(columns=[c for c in merged.columns if c.endswith("_bakta")])

    # replace NaN with "NA" for output
    merged = merged.fillna("NA")

    return merged.to_dict(orient="records")

# ---- summary table ----
def make_summary_table(merged_records):
    df = pd.DataFrame(merged_records)

    def source_type(row):
        has_df = row.get("df_gene_name") not in (None, "NA", "")
        has_pl = row.get("pl_gene_name") not in (None, "NA", "")
        if has_df and has_pl:
            return "Both"
        elif has_df:
            return "DF only"
        elif has_pl:
            return "Padloc only"
        else:
            return "NA"

    summary = pd.DataFrame()
    summary["locus_tag"] = df["locus_tag"]
    summary["source_type"] = df.apply(source_type, axis=1)
    summary["consolidated_gene"] = df.apply(lambda r: r["df_gene_name"] if r["df_gene_name"] not in (None, "NA", "") else r["pl_gene_name"], axis=1)
    summary["consolidated_system"] = df.apply(lambda r: r["df_system"] if r["df_system"] not in (None, "NA", "") else r["pl_system"], axis=1)

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


# ---- main ----
def main():
    args = build_parser().parse_args()

    # load + tidy
    df_rows = tidy_df_genes(load_df_genes(args.defensefinder))
    pl_rows = tidy_padloc(load_padloc(args.padloc))
    bakta_rows = load_bakta_clean(args.bakta)

    merged = merge_all(df_rows, pl_rows, bakta_records=bakta_rows)

    outpath = args.outdir / "df_padloc_merged.tsv"
    # define field order with coordinates last
    if merged:
        cols = [c for c in merged[0].keys() if c not in ("start", "end", "strand")]
        field_order = cols + ["start", "end", "strand"]
    else:
        field_order = ["locus_tag", "df_model", "df_gene_name", "df_system", "df_hit_i_eval", 
                        "df_hit_profile_cov", "df_hit_seq_cov", "df_hit_status", "df_sys_wholeness",
                        "df_hit_score", "sample_name", "pl_system", "pl_gene_name", "pl_evalue",
                        "pl_domain_ievalue", "pl_target_cov", "pl_hmm_cov", "start", "end", "strand"]
    write_tsv(merged, outpath, field_order=field_order)

    summary = make_summary_table(merged)
    summary_path = args.outdir / "df_padloc_consolidated.tsv"
    write_tsv(summary, summary_path, field_order=["locus_tag", "source_type", "consolidated_gene", "consolidated_system"])
    
if __name__ == "__main__":
    main()
