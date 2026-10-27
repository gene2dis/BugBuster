"""ARGs-OAP (SARG) tables of bin/results_summary.py, from ARGS_OAP_REPORT's args_oap_long.tsv and
args_oap_cells.tsv.

ARGS_OAP_REPORT (bin/args_oap_report.py) sums each sample's ARGs-OAP rows (pairs and singletons) and recomputes
stage two's normalisation over them: copies_per_cell = SARG copies / nCell (ARGs-OAP's cell estimate from
essential single-copy genes), copies_per_16S = copies / n16S, ppm = reads per million single reads. A copy is a
read weighted by the share of the reference it covers, so copies per cell is not RPKM and not comparable to the
K-MARVEL / RGI figures. SARG names its types and subtypes differently from CARD: ARGs-OAP is reported in its own
section, not matched gene by gene with the CARD methods, and has no tiers.
"""

import numpy as np
import pandas as pd

from summary_common import NA_NOT_RUN, die, fill_columns

LONG_NEEDED = {"sample", "level", "sarg_type", "sarg_subtype", "sarg_gene", "read_count", "copies",
               "copies_per_cell", "copies_per_16S", "ppm"}
CELLS_NEEDED = {"sample", "nRead", "n16S", "nCell", "status"}
NOT_NORMALISED = "not_normalised"


def _read(path, needed):
    df = pd.read_csv(path, sep="\t", dtype={"sample": str})
    missing = needed - set(df.columns)
    if missing:
        die(f"{path} lacks columns: {', '.join(sorted(missing))}")
    return df


def load_args_oap(root):
    """-> (SARG hits at subtype and type level, per-sample counts, samples ARGs-OAP covered), or Nones."""
    sdir = root / "05_arg_prediction" / "read_level" / "args_oap_summary"
    lpath, cpath = sdir / "args_oap_long.tsv", sdir / "args_oap_cells.tsv"
    if not cpath.exists():
        return None, None, set()
    cells = _read(cpath, CELLS_NEEDED)
    long = _read(lpath, LONG_NEEDED) if lpath.exists() else pd.DataFrame(columns=sorted(LONG_NEEDED))
    long = long[long["level"].isin(["type", "subtype"])].reset_index(drop=True)
    return long, cells, set(cells["sample"])


def sarg_status(status, cells):
    """sample_status() plus not_normalised: ARGs-OAP found no 16S or no cells, so nothing could be normalised."""
    out = dict(status)
    for s in cells.loc[cells["status"] == NOT_NORMALISED, "sample"]:
        if out.get(s) == "ok":
            out[s] = NOT_NORMALISED
    return out


def sarg_tables(long, order, status):
    """Wide copies-per-cell matrices per SARG type and per subtype (with its type). A cell for an entry not found
    is 0, NA_low_depth or NA_not_run; blank for a sample that could not be normalised."""
    def fill(mat):
        fills = {s: 0 if status[s] == "ok" else "" if status[s] == NOT_NORMALISED else status[s] for s in order}
        return fill_columns(mat.reindex(columns=order), fills)

    out = {}
    for level, keys, name in (("type", ["sarg_type"], "arg_sarg_type_cpc"),
                              ("subtype", ["sarg_type", "sarg_subtype"], "arg_sarg_subtype_cpc")):
        rows = long[long["level"] == level]
        if rows.empty:
            continue
        g = rows.groupby(keys + ["sample"])["copies_per_cell"]
        mat = g.sum(min_count=1).unstack("sample")
        # a hit without copies per cell (no cell estimate) is blank: neither 0 nor 'not detected'
        blank = (g.size().unstack("sample").notna() & mat.isna()).reindex(columns=order, fill_value=False)
        out[name] = fill(mat).mask(blank, "").sort_index()
    return out


def sarg_sample_counts(long, cells, order, status):
    """Per sample: ARGs-OAP 16S and cell estimates, distinct SARG subtypes and their summed copies per cell."""
    out = pd.DataFrame(index=order)
    c = cells.drop_duplicates("sample").set_index("sample").reindex(order)
    out["args_oap_n16S"] = c["n16S"]
    out["args_oap_ncell"] = c["nCell"]
    sub = long[long["level"] == "subtype"]
    out["sarg_n_subtypes"] = sub.groupby("sample")["sarg_subtype"].nunique().reindex(order).fillna(0).astype(int)
    out["sarg_copies_per_cell"] = (long[long["level"] == "type"].groupby("sample")["copies_per_cell"]
                                   .sum(min_count=1).reindex(order).fillna(0))
    blank = [s for s in order if status[s] in (NA_NOT_RUN, NOT_NORMALISED)]
    out.loc[blank, ["sarg_n_subtypes", "sarg_copies_per_cell"]] = np.nan
    out.loc[[s for s in order if status[s] == NA_NOT_RUN], ["args_oap_n16S", "args_oap_ncell"]] = np.nan
    return out
