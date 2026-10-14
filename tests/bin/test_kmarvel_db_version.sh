#!/bin/bash
#
# Regression tests for bin/kmarvel_db_version.sh, which gives the KMARVEL
# versions.yml its "card:" value:
#   - DB_VERSION present (FORMAT_KMARVEL_DB builds) -> its first line
#   - DB_VERSION missing, empty or blank (--custom_kmarvel_db) -> "custom"
#     plus the md5 of CARD_protein_all.fasta
#   - characters that would break the double-quoted YAML value are dropped
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
SCRIPT="${REPO_DIR}/bin/kmarvel_db_version.sh"

TMP_DIR=$(mktemp -d)
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

PASS=0
FAIL=0

check() {
    # check <desc> <command...>
    local desc="$1"
    shift
    if "$@"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}"
        FAIL=$((FAIL + 1))
    fi
}

fails() {
    ! "$@" > /dev/null 2>&1
}

cd "${TMP_DIR}"
mkdir db
printf '>ARO:3000873|TEM-1\nMSIQHFRVALIPFFAAFCLPVFA\n' > db/CARD_protein_all.fasta
md5=$(md5sum db/CARD_protein_all.fasta | cut -d' ' -f1)
fallback="custom (CARD_protein_all.fasta md5 ${md5})"

# FORMAT_KMARVEL_DB layout
printf '%s\n' 'CARD 4.0.1 (broadstreet-v4.0.1, https://card.mcmaster.ca)' > db/DB_VERSION
check "DB_VERSION present: printed as is" \
    test "$("${SCRIPT}" db)" = 'CARD 4.0.1 (broadstreet-v4.0.1, https://card.mcmaster.ca)'

printf 'CARD 4.0.1\nsecond line\n' > db/DB_VERSION
check "multi-line DB_VERSION: first line only" test "$("${SCRIPT}" db)" = 'CARD 4.0.1'

printf 'CARD "4.0.1" \\ build\r\n' > db/DB_VERSION
check "quotes, backslash and CR dropped" test "$("${SCRIPT}" db)" = 'CARD 4.0.1  build'

# --custom_kmarvel_db without provenance
rm db/DB_VERSION
check "DB_VERSION missing: custom + FASTA md5" test "$("${SCRIPT}" db)" = "${fallback}"

: > db/DB_VERSION
check "DB_VERSION empty: custom + FASTA md5" test "$("${SCRIPT}" db)" = "${fallback}"

printf '  \n' > db/DB_VERSION
check "DB_VERSION blank: custom + FASTA md5" test "$("${SCRIPT}" db)" = "${fallback}"

check "fallback md5 follows the FASTA content" \
    sh -c 'printf ">x\nMK\n" >> db/CARD_protein_all.fasta; test "$("$1" db)" != "$2"' _ "${SCRIPT}" "${fallback}"

check "symlinked FASTA (staged input): md5 of its target" \
    sh -c 'mkdir linked && ln -s "$PWD/db/CARD_protein_all.fasta" linked/CARD_protein_all.fasta &&
           test "$("$1" linked)" = "custom (CARD_protein_all.fasta md5 $(md5sum db/CARD_protein_all.fasta | cut -d" " -f1))"' _ "${SCRIPT}"

mkdir nofasta
check "no FASTA and no DB_VERSION: exits non-zero" fails "${SCRIPT}" nofasta
check "no argument: exits non-zero" fails "${SCRIPT}"

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
