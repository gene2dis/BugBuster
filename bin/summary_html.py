"""Table and Excel output of bin/results_summary.py (the HTML report is bin/summary_site.py)."""

import pandas as pd

from summary_common import warn
from summary_docs import table_doc


def write_tables(outdir, tables):
    tdir = outdir / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    for name, (df, index) in tables.items():
        df.to_csv(tdir / f"{name}.tsv", sep="\t", index=index)


def write_excel(outdir, tables, group_by=None, units=None):
    engine = None
    for mod in ("openpyxl", "xlsxwriter"):
        try:
            __import__(mod)
            engine = mod
            break
        except ImportError:
            continue
    if engine is None:
        warn("neither openpyxl nor xlsxwriter is installed: no Excel workbook written")
        return None
    path = outdir / "results_summary.xlsx"
    # README sheet: every table written, its description, then its documented columns
    readme = []
    for name in tables:
        doc = table_doc(name, group_by)
        if doc is None:
            continue
        unit = (units or {}).get(name)
        readme.append({"table": name[:31], "column": "", "description": doc["title"] + (f" ({unit})" if unit else "")
                       + ": " + doc["description"]})
        for col, desc in doc["columns"]:
            readme.append({"table": name[:31], "column": col, "description": desc})
    with pd.ExcelWriter(path, engine=engine) as xw:
        pd.DataFrame(readme).to_excel(xw, sheet_name="README", index=False)
        for name, (df, index) in tables.items():
            if len(df) > 1_000_000:
                warn(f"{name} has {len(df)} rows: too many for Excel, see tables/{name}.tsv")
                continue
            out = df.reset_index() if index else df
            out.to_excel(xw, sheet_name=name[:31], index=False, freeze_panes=(1, 0))
    return path
