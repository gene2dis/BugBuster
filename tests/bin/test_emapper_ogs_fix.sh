#!/bin/bash
#
# Tests for bin/emapper_ogs_fix.py (design doc Q23): the runtime patch for
# eggNOG-mapper 3.0.0-beta6's OG-string parser (upstream issue #620).
# Asserts:
#   - eggNOG 7 OG names whose family contains '|' keep their full key and
#     taxid level (the issue's fixtures, with and without a taxon name),
#   - normal OG names and multi-OG strings parse as before,
#   - entries that do not fit the format fall back to the original split,
#   - the version guard refuses to patch a method that is not the beta6 one.
#
# Runs on the host (python3 only). With EMAPPER_RUN set to a container
# prefix (e.g. "apptainer exec eggnog-mapper-3.0.0-beta6.sif" or
# "docker run --rm -v $PWD:$PWD -w $PWD <image>"), it also checks the guard
# accepts the installed beta6 method and that the patched class parses the
# issue's fixtures inside that container.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
FIX="${REPO_DIR}/bin/emapper_ogs_fix.py"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

echo "=== emapper_ogs_fix.py test suite ==="

# Written to a file (not stdin) so inspect.getsource works on the test's own
# classes and the guard's hash-mismatch path is exercised
cat > "${TMP_DIR}/host_tests.py" <<'EOF'
import importlib.util, sys

spec = importlib.util.spec_from_file_location("emapper_ogs_fix", sys.argv[1])
fix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fix)
parse = lambda s: fix.parse_ogs_string(None, s)

passed = failed = 0
def check(desc, got, want):
    global passed, failed
    if got == want:
        print(f"✓ {desc}"); passed += 1
    else:
        print(f"✗ {desc}: got {got!r}, want {want!r}"); failed += 1

# Upstream #620 fixtures
check("pipe in family name keeps the full OG key",
      parse("ABC_tran|TL31Y9@131567|A-1"), [("ABC_tran|TL31Y9@131567|A-1", "131567")])
check("pipe in family name with trailing taxon name",
      parse("ABC_tran|TL31Y9@131567|A-1|cellular organisms"),
      [("ABC_tran|TL31Y9@131567|A-1", "131567")])
check("normal family name (control)",
      parse("ABC_membrane_2@2759|G-3!|Eukaryota"), [("ABC_membrane_2@2759|G-3!", "2759")])
check("multi-OG string as stored in eggnog.db (order kept)",
      parse("ABC_membrane_2@131567|A-1,ABC_membrane_2@2759|G-3!,ABC_tran|TL31Y9@131567|A-1"),
      [("ABC_membrane_2@131567|A-1", "131567"), ("ABC_membrane_2@2759|G-3!", "2759"),
       ("ABC_tran|TL31Y9@131567|A-1", "131567")])
check("blank entries skipped", parse(" , PLDc_2@131567|FK-11 ,"), [("PLDc_2@131567|FK-11", "131567")])
# Fallback: entries that do not fit family@taxid|clade keep the original split
check("no '@': original split, level '-'", parse("famA|cladeB"), [("famA|cladeB", "-")])
check("non-numeric taxid: original split", parse("fam@abc|X"), [("fam@abc|X", "abc")])
check("single field, no clade: dropped as before", parse("fam@131567"), [])

# Version guard: a method that is not the beta6 one must not be patched
class Changed:
    def _parse_ogs_string(self, ogs_string):
        return []
try:
    fix.patch(Changed)
    print("✗ guard refuses a changed method: patched anyway"); failed += 1
except SystemExit as e:
    ok = "not the v3.0.0-beta6 method" in str(e)
    print(("✓" if ok else "✗") + " guard refuses a changed method")
    passed += ok; failed += (not ok)
check("refused method left untouched", Changed()._parse_ogs_string("x@1|A"), [])

# A method whose source cannot be read (no file behind it) is refused too
ns = {}
exec("class NoSource:\n    def _parse_ogs_string(self, s):\n        return []\n", ns)
try:
    fix.patch(ns["NoSource"])
    print("✗ guard refuses a method without source: patched anyway"); failed += 1
except SystemExit:
    print("✓ guard refuses a method without source"); passed += 1

print(f"\n{passed} passed, {failed} failed")
sys.exit(1 if failed else 0)
EOF
python3 -B "${TMP_DIR}/host_tests.py" "${FIX}"

if [ -n "${EMAPPER_RUN:-}" ]; then
    echo
    echo "=== inside the eggNOG-mapper container (${EMAPPER_RUN}) ==="
    # shellcheck disable=SC2086
    ${EMAPPER_RUN} python3 -B - "${FIX}" <<'EOF'
import importlib.util, sys
from eggnogmapper.annotator.e7.annotate import AnnotationEngine

spec = importlib.util.spec_from_file_location("emapper_ogs_fix", sys.argv[1])
fix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fix)

orig = AnnotationEngine._parse_ogs_string(None, "ABC_tran|TL31Y9@131567|A-1")
assert orig == [("ABC_tran|TL31Y9@131567", "-")], orig   # the #620 bug is present
print("✓ installed beta6 method truncates the pipe-family OG (bug present)")
fix.patch(AnnotationEngine)                               # guard accepts beta6
print("✓ guard accepts the installed beta6 method")
got = AnnotationEngine._parse_ogs_string(None, "ABC_tran|TL31Y9@131567|A-1")
assert got == [("ABC_tran|TL31Y9@131567|A-1", "131567")], got
print("✓ patched class parses the pipe-family OG correctly")
EOF
fi
