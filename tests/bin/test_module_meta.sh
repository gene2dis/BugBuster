#!/bin/bash
#
# Checks that every local module documents its interface in meta.yml
# (CONTRIBUTING.md "Module meta.yml") and that the descriptor matches the
# process declarations in main.nf:
#   - every modules/local/<name>/ has a meta.yml that parses as YAML,
#   - `name` matches the directory,
#   - per process: the meta.yml input entries are the process inputs (meta
#     maps plus every val/path variable, in order) and the output entries
#     are exactly the process `emit:` names (plus optional meta entries),
#   - every process in the file has a tag and a label.
# Single-process files use top-level input/output; multi-process files use a
# `processes:` list (name, description, input, output).
#
# Runs on the host (python3 + PyYAML).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

echo "=== module meta.yml / declaration consistency ==="

python3 - "${REPO_DIR}/modules/local" <<'PYEOF'
import re, sys
from pathlib import Path
import yaml

root = Path(sys.argv[1])
problems = []

def block(body, start, ends):
    """Text of a process section (input:/output:) up to the next section."""
    m = re.search(rf'^\s*{start}:\s*$', body, re.M)
    if not m:
        return ''
    rest = body[m.end():]
    e = re.search(rf'^\s*(?:{"|".join(ends)}):', rest, re.M)
    return rest[:e.start()] if e else rest

def input_names(text):
    names = []
    for line in text.splitlines():
        line = line.split('//')[0]
        for kind, arg in re.findall(r'\b(val|path|env)\s*\(\s*([^,)\s]+)', line):
            names.append(arg.strip('"\''))
        m = re.match(r'\s*(?:val|path)\s+([A-Za-z_]\w*)\s*$', line)
        if m:
            names.append(m.group(1))
    return names

def emit_names(text):
    return re.findall(r'emit:\s*([A-Za-z_]\w*)', text)

def entry_names(entries):
    out = []
    for e in entries or []:
        if isinstance(e, dict):
            out.extend(e.keys())
        else:
            out.append(str(e))
    return out

for d in sorted(p for p in root.iterdir() if p.is_dir()):
    nf, meta = d / 'main.nf', d / 'meta.yml'
    if not nf.exists():
        continue
    if not meta.exists():
        problems.append(f'{d.name}: no meta.yml'); continue
    try:
        doc = yaml.safe_load(meta.read_text())
    except yaml.YAMLError as err:
        problems.append(f'{d.name}: meta.yml does not parse: {err}'); continue
    if doc.get('name') != d.name:
        problems.append(f"{d.name}: name is {doc.get('name')!r}")
    src = nf.read_text()
    procs = {}
    for m in re.finditer(r'^process\s+(\w+)\s*\{', src, re.M):
        nxt = re.search(r'^process\s+\w+\s*\{', src[m.end():], re.M)
        procs[m.group(1)] = src[m.end(): m.end() + nxt.start()] if nxt else src[m.end():]
    if len(procs) == 1:
        documented = {next(iter(procs)): doc}
    else:
        documented = {p.get('name'): p for p in (doc.get('processes') or [])}
        missing = set(procs) - set(documented)
        extra = set(documented) - set(procs)
        if missing: problems.append(f'{d.name}: processes not documented: {sorted(missing)}')
        if extra: problems.append(f'{d.name}: documented processes not in main.nf: {sorted(extra)}')
    for pname, body in procs.items():
        if not re.search(r'^\s*tag\s', body, re.M):
            problems.append(f'{d.name}/{pname}: no tag')
        if not re.search(r'^\s*label\s', body, re.M):
            problems.append(f'{d.name}/{pname}: no label')
        pd = documented.get(pname)
        if pd is None:
            continue
        want_in = input_names(block(body, 'input', ['output', 'when', 'script', 'exec', 'stub']))
        got_in = entry_names(pd.get('input'))
        if got_in != want_in:
            problems.append(f'{d.name}/{pname}: input {got_in} != declared {want_in}')
        want_out = emit_names(block(body, 'output', ['when', 'script', 'exec', 'stub']))
        # nf-core style lists the meta map of tuple outputs too; not an emit
        got_out = [o for o in entry_names(pd.get('output')) if o not in ('meta', 'meta_list')]
        if sorted(got_out) != sorted(want_out):
            problems.append(f'{d.name}/{pname}: output {sorted(got_out)} != emits {sorted(want_out)}')

checked = sum(1 for d in root.iterdir() if (d / 'main.nf').exists())
for p in problems:
    print(f'✗ {p}')
print(f'\n{checked} modules checked, {len(problems)} problem(s)')
sys.exit(1 if problems else 0)
PYEOF
