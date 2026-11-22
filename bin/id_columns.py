"""Read a table with its id columns as text (roadmap Phase 20, RV2).

pandas infers column types: a sample id ``001`` is read back as the number 1,
``1e5`` as 100000.0, ``True`` as a boolean, and ``NA`` / ``None`` / ``null`` as
a missing value. Ids (samples, bins, assemblies, database labels) must stay
the text the pipeline wrote, or merges on them fail or silently drop rows.

``read_table(path, ids)`` reads the named columns (and, with ``index=True``, the
first column, which becomes the index) as text with no missing-value words,
while every other column keeps pandas' normal parsing, missing values
included: ``keep_default_na=False`` alone would turn the blank cells of
numeric columns into empty strings.
"""

import pandas as pd

# pandas' own list of the strings read as a missing value by default
try:
    from pandas._libs.parsers import STR_NA_VALUES as _NA_VALUES
except ImportError:  # pragma: no cover - older or newer pandas layouts
    _NA_VALUES = {"", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND",
                  "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null"}
NA_VALUES = sorted(_NA_VALUES)


def text_options(columns, ids):
    """read_csv keywords that read the columns in ``ids`` (those present) as
    text without missing values, and leave the others to pandas."""
    ids = [c for c in ids if c in columns]
    return {
        "dtype": {c: str for c in ids},
        "keep_default_na": False,
        "na_values": {c: NA_VALUES for c in columns if c not in ids},
    }


def read_table(path, ids=(), index=False, **kwargs):
    """pd.read_csv(path, **kwargs) with the columns named in ``ids`` read as
    text; with ``index=True`` the first column is the index, also read as
    text. Missing id columns are ignored. kwargs must not set dtype,
    keep_default_na, na_values or index_col."""
    header_kwargs = {k: v for k, v in kwargs.items() if k in ("sep", "comment", "skiprows", "encoding", "header")}
    columns = list(pd.read_csv(path, nrows=0, **header_kwargs).columns)
    ids = list(ids)
    if index and columns:
        ids.append(columns[0])
    options = text_options(columns, ids)
    df = pd.read_csv(path, **options, **kwargs)
    if index and columns:
        df = df.set_index(columns[0])
    return df
