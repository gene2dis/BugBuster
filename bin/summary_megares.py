"""KARGA / KARGVA (MEGARes) tables of bin/results_summary.py, from ARG_NORM_REPORT's KARGA_norm.csv and
KARGVA_norm.csv.

ARG_NORM_REPORT (bin/arg_norm_report.py) keeps KARGA genes with >= 90% of the gene covered by k-mers and
KARGVA variants with >= 2 SNP k-mers and >= 80% coverage, and joins each sample's ARGs-OAP counts (nRead
single reads, n16S, nCell). Its abundance is a k-mer depth: copies_per_cell = AverageKMerDepth / nCell and
CPM_ARGs = AverageKMerDepth per million single reads (not read counts, so not comparable to the RPM / RPKM of
K-MARVEL and RGI). A sample with no passing hit has one row with blank ARG fields.

MEGARes headers are MEG_id|type|class|mechanism|group (KARGA) and id|mutation|ARO|description|MEG_id|type|
class|mechanism|group|... (KARGVA). Alleles of one MEGARes group carry nearly the same k-mer depth (they share
k-mers), so a group is represented by its deepest allele, never by the sum.
"""

import numpy as np
import pandas as pd

from summary_common import NA_NOT_RUN, die, fill_columns


KARGA_NEEDED = {"sample", "nRead", "n16S", "nCell", "GeneIdx"}
# hit columns: older ARG_NORM_REPORT runs left CPM_ARGs and copies_per_cell out of a file where no gene or
# variant passed, so a file without hits may lack them (they are blank there anyway)
KARGA_VALUES = ["PercentGeneCovered", "AverageKMerDepth", "CPM_ARGs", "copies_per_cell"]
KARGA_FIELDS = ["megares_id", "megares_type", "megares_class", "megares_mechanism", "megares_group",
                "snp_confirmation_required"]
KARGVA_FIELDS = ["kargva_id", "kargva_mutation", "aro_id", "kargva_description", "megares_id", "megares_type",
                 "megares_class", "megares_mechanism", "megares_group"]


def _read_norm(path):
    df = pd.read_csv(path, dtype={"sample": str, "GeneIdx": str})
    missing = KARGA_NEEDED - set(df.columns)
    if missing:
        die(f"{path} lacks columns: {', '.join(sorted(missing))}")
    hits = df["GeneIdx"].fillna("").str.strip() != ""
    absent = [c for c in KARGA_VALUES if c not in df.columns]
    if absent and hits.any():
        die(f"{path} lacks columns: {', '.join(absent)}")
    for c in absent:  # a file with no hit
        df[c] = np.nan
    df["sample"] = df["sample"].astype(str)
    return df


def _megares(header, kargva):
    """MEGARes header -> fields. KARGA: MEG_id|type|class|mechanism|group[|RequiresSNPConfirmation];
    KARGVA: id|mutation|ARO|description|MEG_id|type|class|mechanism|group|..."""
    parts = str(header).lstrip(">").split("|")
    if kargva:
        parts += [""] * (9 - len(parts))
        return {"kargva_id": parts[0], "kargva_mutation": parts[1], "aro_id": parts[2],
                "kargva_description": parts[3], "megares_id": parts[4], "megares_type": parts[5],
                "megares_class": parts[6], "megares_mechanism": parts[7], "megares_group": parts[8]}
    parts += [""] * (5 - len(parts))
    return {"megares_id": parts[0], "megares_type": parts[1], "megares_class": parts[2],
            "megares_mechanism": parts[3], "megares_group": parts[4],
            "snp_confirmation_required": "RequiresSNPConfirmation" in parts[5:]}


def load_karga(root):
    """-> (KARGA hits, KARGVA hits, ARGs-OAP counts per sample, samples ARG_NORM_REPORT covered), or Nones."""
    sdir = root / "05_arg_prediction" / "read_level" / "summary"
    kpath, vpath = sdir / "KARGA_norm.csv", sdir / "KARGVA_norm.csv"
    if not kpath.exists() and not vpath.exists():
        return None, None, None, set()
    frames, cells = {}, []
    for key, path in (("karga", kpath), ("kargva", vpath)):
        if not path.exists():
            frames[key] = None
            continue
        df = _read_norm(path)
        cells.append(df[["sample", "nRead", "n16S", "nCell"]])
        hits = df[df["GeneIdx"].fillna("").str.strip() != ""].copy()
        # columns given, so a file with no hit still has every field (an empty table, not a missing column)
        fields = pd.DataFrame([_megares(h, key == "kargva") for h in hits["GeneIdx"]], index=hits.index,
                              columns=KARGVA_FIELDS if key == "kargva" else KARGA_FIELDS)
        out = pd.concat([hits[["sample"]], fields], axis=1)
        if key == "kargva":
            out["gene"] = hits.get("Gen", "")
            out["kargva_drug"] = hits.get("Drug_resistance", "")
        out["percent_covered"] = pd.to_numeric(hits["PercentGeneCovered"], errors="coerce")
        out["kmer_depth"] = pd.to_numeric(hits["AverageKMerDepth"], errors="coerce")
        out["copies_per_cell"] = pd.to_numeric(hits["copies_per_cell"], errors="coerce")
        out["cpm_kmer"] = pd.to_numeric(hits["CPM_ARGs"], errors="coerce")
        frames[key] = out.reset_index(drop=True)
    cells = pd.concat(cells).drop_duplicates("sample").rename(
        columns={"nRead": "args_oap_nread", "n16S": "args_oap_n16S", "nCell": "args_oap_ncell"})
    return frames["karga"], frames["kargva"], cells, set(cells["sample"])


def karga_representatives(karga):
    """group_representative: the deepest allele of each MEGARes group in a sample (alleles share k-mers)."""
    if karga.empty:
        karga["group_representative"] = pd.Series(dtype=bool)
        return karga
    karga = karga.sort_values(["sample", "megares_group", "kmer_depth"], ascending=[True, True, False]).copy()
    karga["group_representative"] = ~karga.duplicated(["sample", "megares_group"])
    return karga


def karga_tables(karga, order, status):
    """Wide copies-per-cell matrices from the group representatives: per MEGARes group (with type and class)
    and per class (sum of its groups). A cell for a group not found is 0, NA_low_depth or NA_not_run."""
    rep = karga[karga["group_representative"]]

    def fill(mat):
        return fill_columns(mat.reindex(columns=order), {s: 0 if status[s] == "ok" else status[s] for s in order})

    out = {}
    if rep.empty:
        return out

    def wide(keys):
        g = rep.groupby(keys + ["sample"])["copies_per_cell"]
        mat = g.sum(min_count=1).unstack("sample")
        # a hit without copies per cell (ARGs-OAP nCell 0) is blank: neither 0 nor 'not detected'
        blank = (g.size().unstack("sample").notna() & mat.isna()).reindex(columns=order, fill_value=False)
        return fill(mat).mask(blank, "")
    out["arg_karga_group_cpc"] = wide(["megares_type", "megares_class", "megares_group"])
    out["arg_karga_class_cpc"] = wide(["megares_type", "megares_class"])
    return out


def karga_sample_counts(karga, kargva, order, status):
    """Per sample: distinct MEGARes groups, their summed copies per cell, and KARGVA variants found."""
    out = pd.DataFrame(index=order)
    rep = karga[karga["group_representative"]] if karga is not None else pd.DataFrame(columns=["sample"])
    if karga is not None:
        out["karga_n_groups"] = rep.groupby("sample")["megares_group"].nunique().reindex(order).fillna(0).astype(int)
        out["karga_copies_per_cell"] = rep.groupby("sample")["copies_per_cell"].sum().reindex(order).fillna(0)
    if kargva is not None:
        out["kargva_n_variants"] = kargva.groupby("sample").size().reindex(order).fillna(0).astype(int)
    not_run = [s for s in order if status[s] == NA_NOT_RUN]
    out.loc[not_run] = np.nan
    return out
