#!/bin/bash
#
# Regenerate the committed SUPER-FOCUS read-branch fixture (design doc T8b):
# tests/data/superfocus/{reads,sf_db,sftest.superfocus_*,superfocus_versions.yml,
# sftest.ags.tsv}. Reads come from make_superfocus_fixtures.py (deterministic);
# this wrapper builds the DB_90 fixture databases for BOTH aligners from the
# committed source proteins in the pipeline's pinned SUPER-FOCUS image
# (diamond makedb / mmseqs createdb outputs verified byte-stable), then runs
# SUPER-FOCUS + bin/superfocus_function_profile.py exactly as the SUPERFOCUS
# module does, so AGGREGATE_READ_FUNCTIONS / bin-script tests need no aligner.
#
# Database root layout (what `superfocus -b` expects; FORMAT_SUPERFOCUS_DB
# output and the figshare databases unpack the same way):
#   sf_db/db/database_PKs.txt
#   sf_db/db/static/diamond/90_clusters.db.dmnd
#   sf_db/db/static/mmseqs2/90_clusters.db*
#   sf_db/DB_VERSION
#
# Usage: tests/bin/make_superfocus_fixtures.sh   (needs python3 and docker)

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FIX="${REPO}/tests/data/superfocus"
SF_IMG=$(grep -o "community.wave.seqera.io/library/super-focus[^']*" \
    "${REPO}/modules/local/superfocus/main.nf" | grep -v '^oras' | tail -1)

python3 "${REPO}/tests/bin/make_superfocus_fixtures.py"

rm -rf "${FIX}/sf_db/db/static"
mkdir -p "${FIX}/sf_db/db/static/diamond" "${FIX}/sf_db/db/static/mmseqs2"
printf '%s\n' "SUPER-FOCUS test fixture (tests/bin/make_superfocus_fixtures.sh)" > "${FIX}/sf_db/DB_VERSION"

docker run --rm -u "$(id -u):$(id -g)" -v "${FIX}:/f" -v "${REPO}/bin:/b:ro" -w /f "${SF_IMG}" bash -c '
    set -euo pipefail
    diamond makedb --in source/proteins.faa --db sf_db/db/static/diamond/90_clusters.db --quiet
    (cd sf_db/db/static/mmseqs2 && mmseqs createdb ../../../../source/proteins.faa 90_clusters.db > /dev/null)

    work=$(mktemp -d)
    gzip -cdf reads/sftest_R1.fastq.gz reads/sftest_R2.fastq.gz reads/sftest_S.fastq.gz > "${work}/sftest.fastq"
    n_reads=$(( $(wc -l < "${work}/sftest.fastq") / 4 ))
    cd "${work}"
    superfocus -q sftest.fastq -dir out -a diamond -db DB_90 -b /f/sf_db -t 1 -tmp ./tmp -d -l sf.log
    cp out/output_all_levels_and_function.xls /f/sftest.superfocus_all_levels_and_function.xls
    cd /f
    python3 /b/superfocus_function_profile.py --table sftest.superfocus_all_levels_and_function.xls \
        --query sftest.fastq --input-reads "${n_reads}" --aligner diamond --database 90 \
        --sample-id sftest --prefix sftest
    rm -rf "${work}"'

printf '"SUPERFOCUS":\n    superfocus: 1.8\n    diamond: 2.2.1\n    mmseqs2: 18.8cc5c\n    superfocus_aligner: diamond\n    superfocus_db: SUPER-FOCUS test fixture (tests/bin/make_superfocus_fixtures.sh)\n' \
    > "${FIX}/superfocus_versions.yml"
sed 's/^wtest\t/sftest\t/' "${REPO}/tests/data/woltka/wtest.ags.tsv" > "${FIX}/sftest.ags.tsv"

find "${FIX}" -type f | sort | xargs md5sum
