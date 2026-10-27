#!/usr/bin/env python3
"""
ARGs-OAP run-level report (ARGS_OAP_REPORT): the per-sample ARGs-OAP directories of a run as two tables.

Usage: args_oap_report.py [--outdir DIR] <sample>_args_oap [...]

Each <sample>_args_oap directory holds stage one's metadata.txt (nRead single reads, n16S and nCell per
ARGs-OAP sample), ARGS_OAP's stage_two_status.txt (ok, no_sarg_hit or not_normalised) and, when stage two
found SARG hits, blastout.filtered.txt (one row per read and SARG
reference: count, copy = count x aligned length / reference length, and the reference's gene, subtype and
type). ARGs-OAP treats the singleton file as a sample of its own (<id>_tmp_reads_S) and normalises each by its
own counts; here a BugBuster sample is the sum of its rows, so its normalisation recomputes stage two's
formulas over the pairs and the singletons together:
  copies_per_cell = sum(copy) / nCell    (stage two's normalized_cell)
  copies_per_16S  = sum(copy) / n16S     (normalized_16S)
  ppm             = sum(count) / nRead x 1e6   (ppm; per million single reads)
A zero denominator leaves the cell blank.

Outputs (in --outdir):
  args_oap_long.tsv   sample, level (type | subtype | gene), sarg_type, sarg_subtype, sarg_gene, read_count,
                      copies, copies_per_cell, copies_per_16S, ppm; one row per sample and SARG entry found
  args_oap_cells.tsv  sample, nRead, n16S, nCell, sarg_read_count, status: ok, no_sarg_hit (nothing left after
                      the SARG search) or not_normalised (no 16S or no cells, so stage two stopped)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

SUFFIX = "_args_oap"
LEVELS = [("type", ["type"]), ("subtype", ["type", "subtype"]), ("gene", ["type", "subtype", "gene"])]
LONG_COLUMNS = ["sample", "level", "sarg_type", "sarg_subtype", "sarg_gene", "read_count", "copies",
                "copies_per_cell", "copies_per_16S", "ppm"]
CELL_COLUMNS = ["sample", "nRead", "n16S", "nCell", "sarg_read_count", "status"]
HIT_COLUMNS = {"count", "copy", "gene", "subtype", "type"}


def die(msg):
    sys.exit(f"ERROR: {msg}")


def per_sample(path):
    """-> (cells row, hits per level) for one <sample>_args_oap directory."""
    sample = path.name[:-len(SUFFIX)]
    meta_file = path / "metadata.txt"
    if not meta_file.exists():
        die(f"{meta_file} not found")
    meta = pd.read_csv(meta_file, sep="\t", dtype={"sample": str})
    if meta.empty:
        die(f"{meta_file} has no sample row")
    cells = {c: float(pd.to_numeric(meta[c]).sum()) for c in ("nRead", "n16S", "nCell")}
    cells["nRead"] = int(cells["nRead"])
    hits_file = path / "blastout.filtered.txt"
    hits = pd.read_csv(hits_file, sep="\t") if hits_file.exists() else pd.DataFrame(columns=sorted(HIT_COLUMNS))
    missing = HIT_COLUMNS - set(hits.columns)
    if missing:
        die(f"{hits_file} lacks columns: {', '.join(sorted(missing))}")
    status_file = path / "stage_two_status.txt"
    if status_file.exists():
        status = status_file.read_text().strip()
        if status not in ("ok", "no_sarg_hit", "not_normalised"):
            die(f"{status_file}: unknown status '{status}'")
    elif not hits.empty:
        status = "ok"
    else:  # a directory written before ARGS_OAP recorded the outcome
        status = "not_normalised" if 0 in (cells["nRead"], cells["n16S"], cells["nCell"]) else "no_sarg_hit"
    rows = []
    for level, keys in LEVELS:
        if hits.empty:
            break
        agg = hits.groupby(keys, as_index=False)[["count", "copy"]].sum()
        agg.insert(0, "level", level)
        rows.append(agg)
    long = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["level", "count", "copy"])
    long = long.rename(columns={"type": "sarg_type", "subtype": "sarg_subtype", "gene": "sarg_gene",
                                "count": "read_count", "copy": "copies"})
    long.insert(0, "sample", sample)

    def ratio(num, den, scale=1.0):
        return num / den * scale if den else float("nan")
    long["copies_per_cell"] = long["copies"].map(lambda v: ratio(v, cells["nCell"]))
    long["copies_per_16S"] = long["copies"].map(lambda v: ratio(v, cells["n16S"]))
    long["ppm"] = long["read_count"].map(lambda v: ratio(v, cells["nRead"], 1e6))
    cell_row = {"sample": sample, **cells, "sarg_read_count": float(hits["count"].sum()) if not hits.empty else 0.0,
                "status": status}
    return cell_row, long.reindex(columns=LONG_COLUMNS)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dirs", nargs="+", type=Path, help="<sample>_args_oap directories")
    ap.add_argument("--outdir", type=Path, default=Path("."))
    a = ap.parse_args()
    dirs = []
    for d in a.dirs:
        if not d.is_dir() or not d.name.endswith(SUFFIX):
            die(f"{d} is not a <sample>{SUFFIX} directory")
        dirs.append(d)
    names = [d.name[:-len(SUFFIX)] for d in dirs]
    if len(set(names)) != len(names):
        die("the same sample is given twice")
    cells, longs = [], []
    for d in sorted(dirs, key=lambda p: p.name):
        row, long = per_sample(d)
        cells.append(row)
        longs.append(long)
    level_order = {lvl: i for i, (lvl, _) in enumerate(LEVELS)}
    found = [frame for frame in longs if not frame.empty]
    long = pd.concat(found, ignore_index=True) if found else pd.DataFrame(columns=LONG_COLUMNS)
    long = long.assign(_lvl=long["level"].map(level_order)).sort_values(
        ["sample", "_lvl", "sarg_type", "sarg_subtype", "sarg_gene"], na_position="first").drop(columns="_lvl")
    a.outdir.mkdir(parents=True, exist_ok=True)
    long.to_csv(a.outdir / "args_oap_long.tsv", sep="\t", index=False)
    pd.DataFrame(cells, columns=CELL_COLUMNS).to_csv(a.outdir / "args_oap_cells.tsv", sep="\t", index=False)
    print(f"{len(cells)} samples, {int((long['level'] == 'subtype').sum())} sample x subtype rows")


if __name__ == "__main__":
    main()
