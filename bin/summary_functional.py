"""Read-level functional profile (Woltka, SUPER-FOCUS or HUMAnN) and MicrobeCensus tables of
bin/results_summary.py.

Reads AGGREGATE_READ_FUNCTIONS' tables in 07_functional_annotation/summary/ (read_function_abundance.tsv,
read_annotated_fraction.tsv, read_sample_summary.tsv) and MicrobeCensus' <sample>.ags.tsv. The read branch is
reported on its own, never merged with the contig branch (function_*.tsv), as the pipeline does.

Units: abundance_native is the backend's unit (Woltka and SUPER-FOCUS: reads, mates counted separately and
fractional when a read hits several ORFs; HUMAnN: RPK). abundance_cpge is copies per genome equivalent
(RPK / MicrobeCensus genome equivalents; Woltka and HUMAnN only, blank for a sample without MicrobeCensus).
A read's ORF can carry several terms of one ontology, so the terms of an ontology do not add up to the
annotated reads: the tables are not compositional.
"""

import numpy as np
import pandas as pd

from summary_common import NA_LOW_DEPTH, NA_NO_AGS, NA_NOT_RUN, die, warn


ABUNDANCE_COLUMNS = ["sample_id", "source", "backend", "ontology", "accession", "description", "abundance_tpm",
                     "abundance_cpge", "abundance_native", "native_unit"]

# COG functional categories (NCBI COG 2024); Woltka writes the letters without a description
COG_CATEGORIES = {
    "A": "RNA processing and modification", "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion", "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism", "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only", "S": "Function unknown",
    "T": "Signal transduction mechanisms", "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms", "W": "Extracellular structures", "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure", "Z": "Cytoskeleton",
}


def load_read_functional(root):
    """-> {"abundance", "fraction", "summary", "ran"} from 07_functional_annotation/summary/, or None."""
    sdir = root / "07_functional_annotation" / "summary"
    path = sdir / "read_function_abundance.tsv"
    if not path.exists():
        return None
    ab = pd.read_csv(path, sep="\t", dtype={"sample_id": str, "accession": str, "description": str},
                     keep_default_na=False, na_values=[""])
    missing = [c for c in ABUNDANCE_COLUMNS if c not in ab.columns]
    if missing:
        die(f"{path} lacks columns: {', '.join(missing)}")
    for col in ("abundance_cpge", "abundance_native"):
        ab[col] = pd.to_numeric(ab[col], errors="coerce")
    ab["description"] = ab["description"].fillna("")
    fraction = summary = None
    fpath, spath = sdir / "read_annotated_fraction.tsv", sdir / "read_sample_summary.tsv"
    if fpath.exists():
        fraction = pd.read_csv(fpath, sep="\t", dtype={"sample_id": str})
    if spath.exists():
        summary = pd.read_csv(spath, sep="\t", dtype={"sample_id": str, "db_version": str, "tool_version": str})
    ran = set(ab["sample_id"])
    if summary is not None:
        ran |= set(summary["sample_id"])
    return {"abundance": ab, "fraction": fraction, "summary": summary, "ran": ran}


def load_microbecensus(root):
    """MicrobeCensus <sample>.ags.tsv -> DataFrame(sample, mc_avg_genome_size_bp, mc_genome_equivalents)."""
    mdir = root / "07_functional_annotation" / "microbecensus"
    rows = []
    for path in sorted(mdir.glob("*/*.ags.tsv")) if mdir.is_dir() else []:
        df = pd.read_csv(path, sep="\t", dtype={"sample_id": str})
        if not {"sample_id", "average_genome_size_bp", "genome_equivalents"} <= set(df.columns):
            warn(f"{path}: not a MicrobeCensus AGS table, skipped")
            continue
        for r in df.itertuples(index=False):
            rows.append({"sample": r.sample_id, "mc_avg_genome_size_bp": float(r.average_genome_size_bp),
                         "mc_genome_equivalents": float(r.genome_equivalents)})
    return pd.DataFrame(rows, columns=["sample", "mc_avg_genome_size_bp", "mc_genome_equivalents"])


def merge_read_functional(parts):
    """Batches' read functional tables -> one set; batches must use one backend."""
    parts = [p for p in parts if p is not None]
    if not parts:
        return None
    ab = pd.concat([p["abundance"] for p in parts], ignore_index=True)
    backends = sorted(set(ab["backend"]))
    if len(backends) > 1:
        die(f"results directories use different read functional backends: {backends}")
    cat = lambda key: (pd.concat([p[key] for p in parts if p[key] is not None], ignore_index=True)  # noqa: E731
                       if any(p[key] is not None for p in parts) else None)
    return {"abundance": ab, "fraction": cat("fraction"), "summary": cat("summary"),
            "ran": set().union(*[p["ran"] for p in parts]), "backend": backends[0] if backends else ""}


def function_status(order, ran, low):
    return {s: NA_NOT_RUN if s not in ran else (NA_LOW_DEPTH if s in low else "ok") for s in order}


def function_tables(func, order, status):
    """{ontology: (wide table, unit)}: accession, description and one column per sample, in CPGE when the
    backend computes it (a sample without MicrobeCensus is NA_no_ags), else in the native unit.
    A term not found is 0, NA_low_depth or NA_not_run."""
    ab = func["abundance"]
    has_cpge = ab["abundance_cpge"].notna().any()
    value = "abundance_cpge" if has_cpge else "abundance_native"
    unit = "cpge" if has_cpge else (ab["native_unit"].iloc[0] if not ab.empty else "native")
    no_ags = set()
    if has_cpge:
        per_sample = ab.groupby("sample_id")["abundance_cpge"].apply(lambda x: x.notna().any())
        no_ags = set(per_sample.index[~per_sample])
    out = {}
    for ont, sub in ab.groupby("ontology", sort=True):
        desc = sub.groupby("accession")["description"].agg(lambda x: next((d for d in x if d), ""))
        if ont == "cog":
            desc = desc.where(desc != "", desc.index.map(lambda a: COG_CATEGORIES.get(a, "")))
        # sum(min_count=1): a missing CPGE stays missing (pivot_table's sum would make it 0)
        mat = sub.groupby(["accession", "sample_id"])[value].sum(min_count=1).unstack("sample_id")
        mat = mat.reindex(columns=order).astype(object)
        for s in order:
            if status[s] == "ok" and s in no_ags:
                mat[s] = NA_NO_AGS
            else:
                mat[s] = mat[s].where(mat[s].notna(), 0 if status[s] == "ok" else status[s])
        mat.insert(0, "description", desc.reindex(mat.index).fillna(""))
        num = mat[order].apply(pd.to_numeric, errors="coerce")
        mat = mat.loc[num.mean(axis=1).fillna(0).sort_values(ascending=False, kind="stable").index]
        out[ont] = (mat, unit)
    return out


def function_sample_columns(func, mc, order):
    """samples_qc columns: reads assigned to ORFs (as a share of the backend's reads), annotated fraction,
    MicrobeCensus average genome size and genome equivalents."""
    out = pd.DataFrame(index=order)
    if mc is not None and not mc.empty:
        m = mc.drop_duplicates("sample").set_index("sample")
        out["mc_avg_genome_size_bp"] = m["mc_avg_genome_size_bp"].reindex(order)
        out["mc_genome_equivalents"] = m["mc_genome_equivalents"].reindex(order)
    if func is None:
        return out
    if func["summary"] is not None:
        sm = func["summary"].drop_duplicates("sample_id").set_index("sample_id")
        out["func_reads_assigned"] = pd.to_numeric(sm["reads_assigned"], errors="coerce").reindex(order)
        if "mc_genome_equivalents" not in out:
            for src, dst in (("average_genome_size_bp", "mc_avg_genome_size_bp"),
                             ("genome_equivalents", "mc_genome_equivalents")):
                out[dst] = pd.to_numeric(sm[src], errors="coerce").reindex(order)
    if func["fraction"] is not None:
        fr = func["fraction"]
        fr = fr[fr["ontology"] == "any"].drop_duplicates("sample_id").set_index("sample_id")
        out["func_fraction_annotated"] = pd.to_numeric(fr["fraction_annotated"], errors="coerce").reindex(order)
    return out


def function_long(func):
    """read_function_abundance.tsv of every batch, sample renamed, without the contig-only TPM column."""
    ab = func["abundance"].rename(columns={"sample_id": "sample"})
    cog = (ab["ontology"] == "cog") & (ab["description"] == "")
    ab.loc[cog, "description"] = ab.loc[cog, "accession"].map(COG_CATEGORIES).fillna("")
    return ab.drop(columns=[c for c in ("abundance_tpm", "source") if c in ab.columns])


def cog_matrix(tables, order):
    """COG category x sample numeric matrix (CPGE or native) for the figure, letters with their names."""
    if "cog" not in tables:
        return None, ""
    mat, unit = tables["cog"]
    num = mat[order].apply(pd.to_numeric, errors="coerce")
    num = num[num.fillna(0).gt(0).any(axis=1)]
    if num.empty:
        return None, unit
    labels = [f"{a}  {mat.at[a, 'description']}" for a in num.index]
    num.index = labels
    num = num.sort_index()
    return num.replace(0, np.nan), unit
