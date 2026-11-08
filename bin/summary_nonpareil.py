"""Nonpareil (QC stage, --nonpareil) for bin/results_summary.py: each sample's sequencing coverage, diversity and the
read pairs for 95 % coverage, from the NONPAREIL module's tables in 01_quality_control/nonpareil/<sample>/:
- <sample>_nonpareil.tsv: the Nonpareil R package's curve fit (bin/nonpareil_fit.R), one row;
- <sample>_nonpareil_curve.tsv: the curve as the package reads it from the .npo, and the model's prediction.
Nonpareil works on R1 of one sample's clean reads, so its values hold in any batch.
"""

import numpy as np
import pandas as pd

from summary_common import warn

VALUE_COLUMNS = {"coverage": "nonpareil_coverage", "diversity": "nonpareil_diversity",
                 "depth_95pct_pairs": "nonpareil_depth_95pct"}
NEEDED = {"sample", "reads", "has_model", "warning", *VALUE_COLUMNS}
CURVE_NEEDED = {"sample", "kind", "x_pairs", "coverage", "coverage_p25", "coverage_p75"}


def load_nonpareil(root):
    """{'values': DataFrame (one row per sample), 'curves': DataFrame} from one results directory, or None when
    Nonpareil did not run there."""
    ndir = root / "01_quality_control" / "nonpareil"
    vals, curves = [], []
    for path in sorted(ndir.glob("*/*_nonpareil.tsv")) if ndir.is_dir() else []:
        df = pd.read_csv(path, sep="\t", dtype={"sample": str, "warning": str}, keep_default_na=False,
                         na_values=[""])
        if not NEEDED <= set(df.columns):
            warn(f"{path}: not a Nonpareil table (missing {sorted(NEEDED - set(df.columns))}), skipped")
            continue
        vals.append(df)
        cpath = path.with_name(path.name.replace("_nonpareil.tsv", "_nonpareil_curve.tsv"))
        if cpath.exists():
            c = pd.read_csv(cpath, sep="\t", dtype={"sample": str, "kind": str}, keep_default_na=False,
                            na_values=[""])
            if CURVE_NEEDED <= set(c.columns):
                curves.append(c)
            else:
                warn(f"{cpath}: not a Nonpareil curve table, skipped")
    if not vals:
        return None
    values = pd.concat(vals, ignore_index=True)
    for c in ("reads", *VALUE_COLUMNS):
        values[c] = pd.to_numeric(values[c], errors="coerce")
    values["has_model"] = values["has_model"].astype(str).str.lower() == "true"
    values["warning"] = values["warning"].fillna("")
    curve = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame(columns=sorted(CURVE_NEEDED))
    for c in ("x_pairs", "coverage", "coverage_p25", "coverage_p75"):
        curve[c] = pd.to_numeric(curve[c], errors="coerce")
    return {"values": values, "curves": curve}


def merge_nonpareil(parts):
    """The batches' Nonpareil results as one, or None when it ran in none of them."""
    parts = [p for p in parts if p is not None]
    if not parts:
        return None
    return {"values": pd.concat([p["values"] for p in parts], ignore_index=True),
            "curves": pd.concat([p["curves"] for p in parts], ignore_index=True)}


def sample_columns(npar):
    """samples_qc columns, indexed by sample: coverage (0-1), diversity Nd and read pairs for 95 % coverage."""
    v = npar["values"].drop_duplicates("sample").set_index("sample")
    return v[list(VALUE_COLUMNS)].rename(columns=VALUE_COLUMNS)


def low_coverage(samples, threshold):
    """Samples whose Nonpareil coverage is below threshold (samples without a value are not listed)."""
    if "nonpareil_coverage" not in samples:
        return []
    cov = samples["nonpareil_coverage"]
    return [s for s in samples.index if np.isfinite(cov.get(s, np.nan)) and cov[s] < threshold]
