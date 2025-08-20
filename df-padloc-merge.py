#!/usr/bin/env python3


# ---- imports ---- 
import argparse # parse command-line flags
import csv # for reading the input files
from pathlib import Path # nicer file paths than raw strings. 


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
        "hit_score": "df_hit_score",
        "replicon": "sample_name"
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
        "start": "pl_start",
        "end": "pl_end",
        "strand": "pl_strand"
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
        "Strand": "strand",
        "Product": "product",
        "Type": "type",
        "Gene": "gene",
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
    parser = build_parser()
    args = parser.parse_args()

    # load
    bakta_rows = load_bakta_clean(args.bakta)
    df_raw = load_df_genes(args.defensefinder)
    df_clean = tidy_df_genes(df_raw)
    padloc_raw = load_padloc(args.padloc)
    padloc_clean = tidy_padloc(padloc_raw)

    # write
    write_tsv(bakta_rows,  args.outdir / "bakta_clean.tsv",         ["locus_tag","product","type","gene","start","end","strand"])
    write_tsv(df_clean,    args.outdir / "defensefinder_clean.tsv", ["locus_tag","df_gene_name","df_system","df_model","df_hit_i_eval","df_hit_profile_cov","df_hit_seq_cov","df_hit_status","df_sys_wholeness","df_hit_score","sample_name"])
    write_tsv(padloc_clean, args.outdir / "padloc_kept.tsv",         ["locus_tag","pl_gene_name","pl_system","pl_evalue","pl_domain_ievalue","pl_target_cov","pl_hmm_cov","pl_start","pl_end","pl_strand"])

if __name__ == "__main__":
    main()
