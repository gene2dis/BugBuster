#!/bin/bash
#
# Every tool a module reports in pipeline_info/software_versions.yml must have a row in
# assets/report/tool_references.tsv, the table the report's Methods use to name it and cite it
# (roadmap Phase 19, R12). Without a row, Methods lists the raw version key as a tool with no
# reference, as happened to GUNC and CheckV. Checks:
#   - the version keys of every module: the indented `tool: version` lines of the versions.yml
#     heredocs (cat <<-END_VERSIONS ... END_VERSIONS) in modules/local and modules/nf-core
#     main.nf files (script and stub blocks), the `    tool: %s` lines of printf and the
#     `    tool: {...}` lines of Python f.write calls that write versions.yml, the tool names of
#     nf-core topic tuples (val('<tool>'), eval(...)), and the pipeline / nextflow lines
#     SOFTWARE_VERSIONS adds; a module with a versions.yml output and no key found fails;
#   - a key built at run time (it contains '$') must be listed in DYNAMIC below, with the
#     names it can take;
#   - every key has a row; every row's kind is tool, database or helper; every citekey a row
#     names is an entry of assets/report/references.yml;
#   - every tool and database key has a section (the pipeline step CITATIONS.md lists it
#     under), a row with a section is a key some module still writes, and the committed
#     CITATIONS.md is what tests/bin/make_citations.py writes from the two files (roadmap
#     Phase 29, CL6).
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
# printf '"%s":\\n    tool: %s\\n' / printf '    tool: %s\\n' and Python f.write(f'    tool: {...}\\n')
INLINE = re.compile(r"""(?:\\n|['"])\s{4}([A-Za-z][A-Za-z0-9_.-]*): (?:%s|\{)""")

keys = {}  # key -> first file that writes it


def add(key, where):
    keys.setdefault(key, where)


files = sorted((repo / "modules" / "local").glob("*/main.nf")) + sorted((repo / "modules" / "nf-core").rglob("main.nf"))
for f in files:
    rel = str(f.relative_to(repo))
    inside = False
    found_here = set()
    for line in f.read_text().splitlines():
        if not inside and ("printf" in line or "f.write" in line):
            for m in INLINE.finditer(line):
                add(m.group(1), rel)
                found_here.add(m.group(1))
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
                found_here.add(n)
        else:
            add(key, rel)
            found_here.add(key)
    text = f.read_text()
    if re.search(r'path\s*\(?\s*"versions\.yml"', text) and not found_here and not TOPIC.search(text):
        problems.append(f"{rel}: writes versions.yml, but no version key was found in it (heredoc, printf, f.write)")

sv = repo / "subworkflows" / "local" / "software_versions.nf"
for m in re.finditer(r'"\s{4}([A-Za-z0-9_.-]+):\s', sv.read_text()):
    add(m.group(1), str(sv.relative_to(repo)))

# [ process, tool, version ] tuples mixed into ch_versions by a subworkflow (the
# database releases PREPARE_DATABASES records; roadmap Phase 22, review R2-08):
# [ 'PROCESS', 'tool', ... ]
TUPLE = re.compile(r"""\[\s*['"]([A-Z][A-Z0-9_:]*)['"]\s*,\s*['"]([A-Za-z0-9_.-]+)['"]\s*,""")
tuple_keys = 0
for f in sorted((repo / "subworkflows" / "local").glob("*.nf")):
    for m in TUPLE.finditer(f.read_text()):
        add(m.group(2), str(f.relative_to(repo)))
        tuple_keys += 1
if tuple_keys < 4:
    problems.append(f"only {tuple_keys} versions tuples found in subworkflows/local (expected the four database "
                    "releases of prepare_databases.nf): the TUPLE pattern no longer matches")

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

# CITATIONS.md lists the rows with a section (roadmap Phase 29, CL6): every tool and database a
# module reports has one, and a row with one is still reported (a tool no module runs any more
# keeps its row, for the summary of older results, but leaves CITATIONS.md)
for tool, row in sorted(rows.items()):
    section = (row.get("section") or "").strip()
    if row["kind"] != "helper" and tool in keys and not section:
        problems.append(f"tool_references.tsv: {tool} ({row['kind']}) has no section, so CITATIONS.md leaves it out")
    if section and tool not in keys:
        problems.append(f"tool_references.tsv: {tool} has a section, but no module reports it: clear its section")

print(f"{len(keys)} version keys, {len(rows)} tool_references.tsv rows, {len(refs)} references")
if problems:
    print("\n".join(f"  FAIL: {p}" for p in problems))
    sys.exit(1)
print("  ok: every version key has a tool reference row, and every citekey a references.yml entry")
PYEOF

echo "=== CITATIONS.md vs assets/report/tool_references.tsv + references.yml ==="
python3 "${REPO_DIR}/tests/bin/make_citations.py" --check
