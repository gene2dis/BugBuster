#!/bin/bash
#
# Checks nextflow_schema.json against the config it describes (roadmap
# Phase 25, RV27; review finding R8-12). CI's lint job used to parse the
# schema as JSON only, and skipped even that when the file was missing.
#   1. The schema is a valid JSON Schema (draft 2020-12, the dialect it
#      declares and nf-schema 2.x reads), and every allOf $ref resolves.
#   2. The schema and nextflow.config hold the same parameters: every
#      top-level param of `nextflow config -flat` has a schema entry, every
#      schema entry has a param, and every map-valued param (the database
#      registry maps, read from config/databases.config) is in
#      validation.defaultIgnoreParams (a map missing there stops every run
#      with an "unrecognised parameter" error). A param in that list needs no
#      schema entry (the deprecated bowtie_ref_* registry names, null by default).
#   3. Every default agrees: a schema default equals the config's, and the
#      config value validates against its schema entry (type, enum, pattern).
#      A param without a schema default has a null or empty config value,
#      apart from the COMPUTED ones below.
#
# Requires nextflow and python3 with jsonschema (CI: apt python3-jsonschema).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
SCHEMA="${SCHEMA:-${REPO_DIR}/nextflow_schema.json}"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "=== nextflow_schema.json vs config test suite ==="
(cd "${REPO_DIR}" && nextflow config -flat . > "${TMP_DIR}/flat.txt" 2> "${TMP_DIR}/flat.err") || {
    echo "✗ nextflow config -flat failed:"; tail -5 "${TMP_DIR}/flat.err"; exit 1
}

python3 -I - "${SCHEMA}" "${TMP_DIR}/flat.txt" <<'EOF'
import ast, json, re, sys

from jsonschema import Draft202012Validator
from jsonschema.exceptions import SchemaError, ValidationError

schema_path, flat_path = sys.argv[1], sys.argv[2]

# Params whose default nextflow.config computes at run time, with no schema
# default by design: the two asset paths are built from projectDir, the
# trace suffix from the launch time
COMPUTED = {
    "arg_priority_table": r".*/assets/arg_priority_tiers\.tsv",
    "summary_watchlist": r".*/assets/summary_watchlist\.tsv",
    "trace_report_suffix": r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}",
}

passed = failed = 0
def check(ok, desc, detail=""):
    global passed, failed
    if ok:
        print(f"✓ {desc}"); passed += 1
    else:
        print(f"✗ {desc}" + (f": {detail}" if detail else "")); failed += 1

def value(text):
    # a value as `nextflow config -flat` prints it
    text = text.strip()
    if text == "null":
        return None
    if text in ("true", "false"):
        return text == "true"
    if re.fullmatch(r"-?\d+", text):
        return int(text)
    if re.fullmatch(r"-?\d+\.\d+([eE]-?\d+)?", text):
        return float(text)
    if len(text) >= 2 and text[0] == text[-1] == "'":
        return text[1:-1].replace("\\'", "'").replace("\\\\", "\\")
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return text

schema = json.load(open(schema_path))
try:
    Draft202012Validator.check_schema(schema)
    check(True, "schema is valid JSON Schema (draft 2020-12)")
except SchemaError as e:
    check(False, "schema is valid JSON Schema (draft 2020-12)", e.message)
check(schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema",
      "schema declares draft 2020-12", schema.get("$schema"))

defs = schema.get("$defs", {})
refs = [item.get("$ref", "") for item in schema.get("allOf", [])]
missing = [r for r in refs if not r.startswith("#/$defs/") or r[len("#/$defs/"):] not in defs]
check(not missing, "every allOf $ref names a $defs group", ", ".join(missing))
unused = sorted(set(defs) - {r[len("#/$defs/"):] for r in refs})
check(not unused, "every $defs group is referenced from allOf (its params are validated)", ", ".join(unused))

# Every param is defined in one group only. A param in two groups is
# validated against both entries, and copies drift apart (config_profile_name
# / _description were in generic_options and cloud_hpc_profile_options until
# roadmap Phase 28)
props, dupes = {}, []
entries = [p for g in defs.values() for p in g.get("properties", {}).items()]
entries += list(schema.get("properties", {}).items())
for name, prop in entries:
    if name in props:
        dupes.append(name)
    props[name] = prop
check(not dupes, "no param is defined in two groups", ", ".join(dupes))

flat, nested, ignore = {}, set(), []
for line in open(flat_path):
    key, sep, text = line.rstrip("\n").partition(" = ")
    if not sep:
        continue
    if key == "validation.defaultIgnoreParams":
        ignore = value(text)
    if not key.startswith("params."):
        continue
    key = key[len("params."):]
    if "." in key:
        nested.add(key.split(".")[0])
    else:
        flat[key] = value(text)

check(bool(flat), f"config has top-level params ({len(flat)})")
# a param in validation.defaultIgnoreParams needs no entry: the deprecated
# registry names (bowtie_ref_*, null unless a -c file sets them) are such params
check(sorted(set(flat) - set(props) - set(ignore)) == [], "every config param has a schema entry (or is in validation.defaultIgnoreParams)",
      ", ".join(sorted(set(flat) - set(props) - set(ignore))))
check(sorted(set(props) - set(flat)) == [], "every schema entry is a config param",
      ", ".join(sorted(set(props) - set(flat))))
check(isinstance(ignore, list) and bool(ignore), "validation.defaultIgnoreParams read from the config", repr(ignore))
check(sorted(nested - set(ignore)) == [], "every map-valued param is in validation.defaultIgnoreParams",
      ", ".join(sorted(nested - set(ignore))))
check(sorted(nested & set(props)) == [], "no map-valued param has a schema entry",
      ", ".join(sorted(nested & set(props))))

bad_default, bad_type, no_default = [], [], []
for name in sorted(set(flat) & set(props)):
    prop, cfg = props[name], flat[name]
    if name in COMPUTED:
        if "default" in prop or not (isinstance(cfg, str) and re.fullmatch(COMPUTED[name], cfg)):
            bad_default.append(f"{name} (computed: config {cfg!r}, schema {prop.get('default', '<none>')!r})")
        continue
    if "default" in prop:
        if prop["default"] != cfg or (type(prop["default"]) is bool) != (type(cfg) is bool):
            bad_default.append(f"{name} (config {cfg!r}, schema {prop['default']!r})")
    elif cfg not in (None, ""):
        no_default.append(f"{name} (config {cfg!r})")
    if cfg is not None:
        try:
            Draft202012Validator(prop).validate(cfg)
        except ValidationError as e:
            bad_type.append(f"{name}: {e.message}")

check(not bad_default, "schema defaults equal the config defaults", "; ".join(bad_default))
check(not no_default, "a param without a schema default has no config default either", "; ".join(no_default))
check(not bad_type, "every config default validates against its schema entry", "; ".join(bad_type))

print(f"\n=== Results: {passed} passed, {failed} failed ===")
sys.exit(1 if failed else 0)
EOF
