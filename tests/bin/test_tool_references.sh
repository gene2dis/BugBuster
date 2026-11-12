#!/bin/bash
#
# Every tool a module reports in pipeline_info/software_versions.yml must have a row in
# assets/report/tool_references.tsv, the table the report's Methods use to name it and cite it
# (roadmap Phase 19, R12). Without a row, Methods lists the raw version key as a tool with no
# reference, as happened to GUNC and CheckV. Checks:
#   - the version keys of every module: the indented `tool: version` lines of the versions.yml
#     heredocs (cat <<-END_VERSIONS ... END_VERSIONS) in modules/local and modules/nf-core
#     main.nf files (script and stub blocks), the tool names of nf-core topic tuples
#     (val('<tool>'), eval(...)), and the pipeline / nextflow lines SOFTWARE_VERSIONS adds;
#   - a key built at run time (it contains '$') must be listed in DYNAMIC below, with the
#     names it can take;
#   - every key has a row; every row's kind is tool, database or helper; every citekey a row
#     names is an entry of assets/report/references.yml.
#
# Runs on the host (python3, standard library only).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

echo "=== software_versions.yml keys vs assets/report/tool_references.tsv ==="

python3 - "${REPO_DIR}" <<'PYEOF'
import csv
import re
import sys
from pathlib import Path

repo = Path(sys.argv[1])
problems = []

# keys written at run time: {(module file relative to the repo, key as written): names it can take}
DYNAMIC = {
    ("modules/local/rgi_bwt/main.nf", "${aligner}"): {"kma", "bowtie2", "bwa"},  # --rgi_aligner
}

HEREDOC_START = re.compile(r"cat\s+<<-?\s*END_VERSIONS")
KEY_LINE = re.compile(r"^\s+([^\s:\"][^:\"]*?):(\s|$)")
TOPIC = re.compile(r"""val\(\s*['"]([^'"]+)['"]\s*\)\s*,\s*eval\(""")

keys = {}  # key -> first file that writes it


def add(key, where):
    keys.setdefault(key, where)


files = sorted((repo / "modules" / "local").glob("*/main.nf")) + sorted((repo / "modules" / "nf-core").rglob("main.nf"))
for f in files:
    rel = str(f.relative_to(repo))
    inside = False
    for line in f.read_text().splitlines():
        if not inside:
            if HEREDOC_START.search(line):
                inside = True
            for m in TOPIC.finditer(line):
                add(m.group(1), rel)
            continue
        if line.strip() == "END_VERSIONS":
            inside = False
            continue
        m = KEY_LINE.match(line)
        if not m:
            continue
        key = m.group(1).strip()
        if "$" in key:
            names = DYNAMIC.get((rel, key))
            if names is None:
                problems.append(f"{rel}: version key '{key}' is built at run time; list the names it can take in "
                                "DYNAMIC (tests/bin/test_tool_references.sh)")
                continue
            for n in sorted(names):
                add(n, rel)
        else:
            add(key, rel)

sv = repo / "subworkflows" / "local" / "software_versions.nf"
for m in re.finditer(r'"\s{4}([A-Za-z0-9_.-]+):\s', sv.read_text()):
    add(m.group(1), str(sv.relative_to(repo)))

if len(keys) < 50:
    problems.append(f"only {len(keys)} version keys found: the parsing above no longer matches the modules")

refs = set(re.findall(r"^([A-Za-z0-9_-]+):\s*$", (repo / "assets/report/references.yml").read_text(), re.M))
rows = {}
with open(repo / "assets/report/tool_references.tsv") as fh:
    for row in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
        rows[row["tool"]] = row
        if row["kind"] not in ("tool", "database", "helper"):
            problems.append(f"tool_references.tsv: {row['tool']} has kind '{row['kind']}' (tool, database or helper)")
        for k in filter(None, (row.get("references") or "").split(";")):
            if k not in refs:
                problems.append(f"tool_references.tsv: {row['tool']} cites '{k}', which references.yml does not have")

for key in sorted(keys):
    if key not in rows:
        problems.append(f"version key '{key}' ({keys[key]}) has no row in assets/report/tool_references.tsv")

print(f"{len(keys)} version keys, {len(rows)} tool_references.tsv rows, {len(refs)} references")
if problems:
    print("\n".join(f"  FAIL: {p}" for p in problems))
    sys.exit(1)
print("  ok: every version key has a tool reference row, and every citekey a references.yml entry")
PYEOF
