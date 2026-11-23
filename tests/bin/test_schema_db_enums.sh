#!/bin/bash
#
# Every database-key parameter (--kraken2_db, --bakta_db, ...) has a schema
# enum equal to the keys of its registry map in config/databases.config, so a
# typo such as --bakta_db v6.0-lite is rejected at launch with the valid keys
# listed, instead of a "Cannot get property 'url' on null object" error from
# PREPARE_DATABASES (roadmap Phase 22, RV16; review R1-03: eleven of the
# nineteen had no enum).
#
# The registry keys come from `nextflow config -flat` (the maps nest, e.g.
# SUPER-FOCUS's per-aligner blocks, so the config file is not parsed by hand;
# -flat rather than -o json, which Nextflow 26.04.4 cannot print for this
# config: "No such variable: meta").
#
# Usage: test_schema_db_enums.sh            check, exit 1 on any difference
#        test_schema_db_enums.sh --write    add the missing enums to
#                                           nextflow_schema.json from the registry
#                                           (after the "default" line; an enum
#                                           that differs is reported, not changed)
#
# Requirements: nextflow, python3.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
MODE="${1:-check}"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

(cd "${TMP_DIR}" && nextflow config -flat "${REPO_DIR}" > "${TMP_DIR}/config.flat")

python3 - "${TMP_DIR}/config.flat" "${REPO_DIR}/nextflow_schema.json" "${MODE}" <<'EOF'
import json
import re
import sys

config_path, schema_path, mode = sys.argv[1:4]
# params.<map>.<key>.<field> = ... (a key with a dot or a dash is quoted) and
# params.<name> = '<value>' for the parameters themselves, in file order
params = {}
for line in open(config_path):
    m = re.match(r"^params\.(\w+)\.(?:'([^']+)'|(\w+))\.", line)
    if m:
        params.setdefault(m.group(1), {}).setdefault(m.group(2) or m.group(3), True)
        continue
    m = re.match(r"^params\.(\w+) = '([^']*)'$", line.rstrip("\n"))
    if m:
        params[m.group(1)] = m.group(2)

# database-key parameter -> its registry map in config/databases.config
PAIRS = {
    "phiX_index": "bowtie_ref_genomes_for_build",
    "host_db": "bowtie_ref_host_index",
    "kraken2_db": "kraken_ref_db",
    "sourmash_db": "sourmash_ref_db",
    "sylph_db": "sylph_ref_db",
    "kmarvel_db": "kmarvel_ref_db",
    "metabuli_db": "metabuli_ref_db",
    "genomad_db": "genomad_ref_db",
    "checkv_db": "checkv_ref_db",
    "checkm2_db": "checkm2_ref_db",
    "gunc_db": "gunc_ref_db",
    "gtdbtk_db": "gtdbtk_ref_db",
    "eggnog_db": "eggnog_ref_db",
    "dbcan_db": "dbcan_ref_db",
    "bakta_db": "bakta_ref_db",
    "cog_db": "cog_ref_db",
    "woltka_db": "woltka_ref_db",
    "superfocus_db": "superfocus_ref_db",
    "humann_db": "humann_ref_db",
}

schema = json.load(open(schema_path))
props = {}
for group in schema["$defs"].values():
    props.update(group.get("properties", {}))

failures = []
missing = {}
registry_maps = sorted(k for k in params if k.endswith("_ref_db") or k.startswith("bowtie_ref_"))
unpaired = sorted(set(registry_maps) - set(PAIRS.values()))
if unpaired:
    failures.append(f"registry maps with no parameter in this test: {', '.join(unpaired)}")

for param, regmap in PAIRS.items():
    if regmap not in params:
        failures.append(f"{param}: registry map {regmap} not in the resolved config")
        continue
    keys = list(params[regmap].keys())
    if param not in props:
        failures.append(f"{param}: not in nextflow_schema.json")
        continue
    enum = props[param].get("enum")
    if mode == "--write" and enum is None:
        missing[param] = keys
        continue
    if enum is None:
        failures.append(f"{param}: no enum (registry keys: {', '.join(keys)})")
    elif sorted(enum) != sorted(keys):
        failures.append(f"{param}: enum {enum} != registry keys {keys} ({regmap})")
    elif params.get(param) not in keys:
        failures.append(f"{param}: default {params.get(param)!r} is not a registry key")
    else:
        print(f"✓ {param}: {', '.join(keys)}")

if mode == "--write":
    # Text insertion, so the rest of the file keeps its layout: the enum goes
    # after the property's "default" line, in the multi-line form the schema's
    # other enums use
    lines = open(schema_path).read().split("\n")
    for param, keys in missing.items():
        start = next(i for i, l in enumerate(lines) if l.strip() == f'"{param}": {{')
        dflt = next(i for i in range(start, len(lines)) if lines[i].strip().startswith('"default":'))
        indent = lines[dflt][: len(lines[dflt]) - len(lines[dflt].lstrip())]
        block = [f'{indent}"enum": ['] + [
            f'{indent}    {json.dumps(k)}{"," if n < len(keys) - 1 else ""}' for n, k in enumerate(keys)
        ] + [f"{indent}],"]
        lines[dflt + 1:dflt + 1] = block
    text = "\n".join(lines)
    json.loads(text)
    open(schema_path, "w").write(text)
    print(f"added the enums of {len(missing)} parameters to {schema_path}")
    for f in failures:
        print(f"✗ {f}")
    sys.exit(1 if failures else 0)

for f in failures:
    print(f"✗ {f}")
print(f"\n{len(PAIRS) - len([f for f in failures if not f.startswith('registry maps')])} of {len(PAIRS)} parameters match their registry")
sys.exit(1 if failures else 0)
EOF
