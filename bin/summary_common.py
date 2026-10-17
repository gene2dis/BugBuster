"""Shared constants and helpers of bin/results_summary.py."""

import re
import sys


NA_LOW_DEPTH = "NA_low_depth"
NA_NOT_RUN = "NA_not_run"
NA_NO_AGS = "NA_no_ags"
GTDB_SUFFIX = re.compile(r"_[A-Z]+$")
ACCESSION = re.compile(r"(GC[AF]_\d+\.\d+)")

VERSION_PROCESS = re.compile(r'^\s*"([^"]+)":\s*$')
VERSION_TOOL = re.compile(r"^\s+([A-Za-z0-9_.+\-]+):\s*(\S.*)?$")


def parse_versions(text):
    """software_versions.yml text -> {process: {tool: version}}, line by line.

    Not YAML parsing on purpose: published files are not always valid YAML (KARGVA's tab-indented
    heredoc leaves its terminator in the file; RGI_BWT/RGI_KMER write multi-line junk). A process is a
    quoted key ending in ':' ("QC:FASTP": -> FASTP, the name after the last ':'); its tools are the
    indented 'tool: version' lines below it. A tool with no version, or an error message for one,
    is 'unknown'.
    """
    out, current = {}, None
    for line in str(text or "").splitlines():
        m = VERSION_PROCESS.match(line)
        if m:
            current = out.setdefault(m.group(1).split(":")[-1], {})
            continue
        m = VERSION_TOOL.match(line)
        if m and current is not None:
            value = (m.group(2) or "").strip().strip('"')
            if not value or "invalid option" in value.lower() or "error" in value.lower():
                value = "unknown"
            current.setdefault(m.group(1), value)
    return out


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
