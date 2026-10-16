"""Shared constants and helpers of bin/results_summary.py."""

import re
import sys


NA_LOW_DEPTH = "NA_low_depth"
NA_NOT_RUN = "NA_not_run"
GTDB_SUFFIX = re.compile(r"_[A-Z]+$")
ACCESSION = re.compile(r"(GC[AF]_\d+\.\d+)")


def die(msg):
    sys.exit(f"ERROR: {msg}")


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


def natural_key(value):
    """Sort key that orders V2 before V10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(value))]


def strip_gtdb(name):
    """Remove GTDB placeholder suffixes for display and matching: Enterococcus_B faecium -> Enterococcus faecium."""
    if not isinstance(name, str):
        return name
    return " ".join(GTDB_SUFFIX.sub("", word) for word in name.split(" "))


def norm_taxon(name):
    return strip_gtdb(name).strip().lower() if isinstance(name, str) else ""
