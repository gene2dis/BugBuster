"""Assembly tables of bin/results_summary.py: contig statistics per assembly (BBMap stats.sh on the MEGAHIT
contigs, BBMAP's length filter from CONTIG_FILTER_SUMMARY) and the share of each sample's reads on the contigs
and on their predicted genes (FEATURECOUNTS_GENES' summaries, --contig_level_functional).

Published inputs: 03_assembly/{per_sample/<sample>,coassembly}/<id>_contig.stats,
pipeline_info/contig_filtering_summary.txt, 07_functional_annotation/gene_abundance/<sample>/
<sample>.featureCounts.txt.summary. An assembly is a sample name, or coassembly.
"""

import re

import pandas as pd

from summary_common import NA_NOT_RUN, die

UNITS = {"": 1, "b": 1, "bp": 1, "kb": 1e3, "kbp": 1e3, "mb": 1e6, "mbp": 1e6, "gb": 1e9, "gbp": 1e9}


def _length_bp(text, path):
    """'8.588 Kbp', '129.856 Mb', '465' -> bp."""
    m = re.fullmatch(r"\s*([\d.]+)\s*([A-Za-z]*)\s*", text)
    if not m or m.group(2).lower() not in UNITS:
        die(f"{path}: unexpected length '{text}' (BBMap stats.sh layout changed?)")
    return float(m.group(1)) * UNITS[m.group(2).lower()]


def parse_contig_stats(path):
    """BBMap stats.sh output -> {contigs, length_bp, n50_bp, l50, max_contig_bp, gc}. BBMap writes N/L50 as
    '<count>/<length>' (its 'N50' is the number of contigs, its 'L50' their length), the reverse of the common use:
    the count is l50 here and the length n50."""
    lines = path.read_text().splitlines()
    fields = {}
    for line in lines:
        if ":" in line and line.startswith(("Main genome", "Max contig")):
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip().split("\t")[0].strip()
    if "Main genome contig total" not in fields:
        die(f"{path}: no 'Main genome contig total' line (not a BBMap stats.sh report?)")
    out = {"contigs": int(fields["Main genome contig total"])}
    gc = None
    for i, line in enumerate(lines[:-1]):
        head = line.split("\t")
        if "GC" in head and head[0] == "A":
            gc = float(lines[i + 1].split("\t")[head.index("GC")])
            break
    out["gc"] = gc
    out["length_bp"] = (_length_bp(fields["Main genome contig sequence total"], path)
                        if "Main genome contig sequence total" in fields else None)
    nl = fields.get("Main genome contig N/L50", "")
    if "/" in nl:
        count, length = nl.split("/", 1)
        out["l50"] = int(count)
        out["n50_bp"] = _length_bp(length, path)
    else:
        out["l50"] = out["n50_bp"] = None
    out["max_contig_bp"] = (_length_bp(fields["Max contig length"], path)
                            if "Max contig length" in fields else None)
    return out


def parse_filter_summary(path):
    """CONTIG_FILTER_SUMMARY's 'Sample_ID Contigs_Before Contigs_After Min_Length All_Filtered' block ->
    DataFrame(assembly, contigs_kept, min_contig_length_bp)."""
    rows, inside = [], False
    for line in path.read_text().splitlines():
        if line.startswith("Sample_ID\t"):
            head = line.split("\t")
            inside = True
            continue
        if inside:
            if line.startswith("---"):
                continue
            if not line.strip():
                break
            rows.append(dict(zip(head, line.split("\t"))))
    if not rows:
        return pd.DataFrame(columns=["assembly", "contigs_kept", "min_contig_length_bp"])
    df = pd.DataFrame(rows)
    need = {"Sample_ID", "Contigs_After", "Min_Length"}
    if not need <= set(df.columns):
        die(f"{path}: the filtering table lacks {sorted(need - set(df.columns))}")
    return pd.DataFrame({"assembly": df["Sample_ID"],
                         "contigs_kept": pd.to_numeric(df["Contigs_After"], errors="coerce"),
                         "min_contig_length_bp": pd.to_numeric(df["Min_Length"], errors="coerce")})


def load_assembly(root):
    """-> {"stats": {assembly: parsed stats}, "filter": DataFrame, "counts": {sample: featureCounts summary}} or
    None when no assembly ran."""
    adir = root / "03_assembly"
    paths = sorted(adir.glob("coassembly/*_contig.stats")) + sorted(adir.glob("per_sample/*/*_contig.stats"))
    if not paths:
        return None
    stats = {p.name[:-len("_contig.stats")]: parse_contig_stats(p) for p in paths}
    fpath = root / "pipeline_info" / "contig_filtering_summary.txt"
    filt = parse_filter_summary(fpath) if fpath.exists() else None
    counts = {}
    gdir = root / "07_functional_annotation" / "gene_abundance"
    for p in sorted(gdir.glob("*/*.featureCounts.txt.summary")) if gdir.is_dir() else []:
        df = pd.read_csv(p, sep="\t")
        if df.columns[0] != "Status" or df.shape[1] != 2:
            die(f"{p}: not a featureCounts summary (expected Status and one count column)")
        counts[p.name[:-len(".featureCounts.txt.summary")]] = dict(zip(df["Status"], df.iloc[:, 1]))
    return {"stats": stats, "filter": filt, "counts": counts}


def assembly_stats_table(parts):
    """assembly_stats: one row per (batch, assembly). parts: [(batch label, load_assembly result)]."""
    rows = []
    for label, a in parts:
        if a is None:
            continue
        filt = a["filter"].set_index("assembly") if a["filter"] is not None else None
        for asm, st in a["stats"].items():
            kept = filt.at[asm, "contigs_kept"] if filt is not None and asm in filt.index else None
            minlen = filt.at[asm, "min_contig_length_bp"] if filt is not None and asm in filt.index else None
            rows.append({
                "assembly": asm, "batch": label,
                "assembly_mode": "coassembly" if asm == "coassembly" else "assembly",
                "contigs": st["contigs"],
                "assembly_length_mb": None if st["length_bp"] is None else round(st["length_bp"] / 1e6, 3),
                "n50_kb": None if st["n50_bp"] is None else round(st["n50_bp"] / 1e3, 3),
                "l50": st["l50"],
                "max_contig_kb": None if st["max_contig_bp"] is None else round(st["max_contig_bp"] / 1e3, 3),
                "gc_pct": None if st["gc"] is None else round(st["gc"] * 100, 2),
                "min_contig_length_bp": minlen,
                "contigs_kept": kept,
                "pct_contigs_kept": (round(kept / st["contigs"] * 100, 2)
                                     if kept is not None and st["contigs"] else None),
            })
    if not rows:
        return None
    df = pd.DataFrame(rows)
    for c in ("contigs", "l50", "min_contig_length_bp", "contigs_kept"):
        df[c] = df[c].astype("Int64")
    return df.sort_values(["batch", "assembly"], kind="stable").reset_index(drop=True)


def read_use_frame(parts, batch_samples):
    """Per sample: the assembly its reads were counted against and the % of reads featureCounts placed on the
    contigs and on predicted genes; NA_not_run for every sample of an assembly batch without featureCounts
    summaries (no --contig_level_functional). parts: [(batch label, load_assembly result)];
    batch_samples: {batch label: [samples]}. -> DataFrame(sample, assembly, pct_reads_on_contigs,
    pct_reads_on_genes) or None."""
    rows = []
    for label, a in parts:
        if a is None:
            continue
        coassembly = "coassembly" in a["stats"]
        for s in batch_samples.get(label, []):
            asm = "coassembly" if coassembly else s
            if asm not in a["stats"]:
                continue
            c = a["counts"].get(s)
            if c is None:
                rows.append({"sample": s, "assembly": asm, "pct_reads_on_contigs": NA_NOT_RUN,
                             "pct_reads_on_genes": NA_NOT_RUN})
                continue
            total = sum(c.values())
            rows.append({"sample": s, "assembly": asm,
                         "pct_reads_on_contigs": round((1 - c.get("Unassigned_Unmapped", 0) / total) * 100, 2)
                         if total else None,
                         "pct_reads_on_genes": round(c.get("Assigned", 0) / total * 100, 2) if total else None})
    return pd.DataFrame(rows) if rows else None
