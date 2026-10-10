#!/bin/bash
#
# Regenerate the committed Woltka read-branch fixture (design doc T8a):
# tests/data/woltka/{db,reads,cog.def.tab,wtest.wol.sam.gz,wtest.woltka_*.tsv,
# woltka_versions.yml}. The sequences, coordinates and maps come from
# make_woltka_fixtures.py (deterministic); this wrapper then builds the
# Bowtie2 index in the pipeline's pinned bowtie2 container with
# --large-index, so the files carry the same .bt2l suffix as the real WoLr2
# index (databases/bowtie2/WoLr2.*.bt2l). --ftabchars 4 keeps the lookup
# table tiny (the 10-char default alone is 8 MB per .1/.rev.1 file).
#
# Usage: tests/bin/make_woltka_fixtures.sh   (needs python3 and docker)

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BOWTIE2_IMG="quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0"

python3 "${REPO}/tests/bin/make_woltka_fixtures.py"

DB="${REPO}/tests/data/woltka/db"
rm -rf "${DB}/databases/bowtie2"
mkdir -p "${DB}/databases/bowtie2"
docker run --rm -u "$(id -u):$(id -g)" -v "${DB}:/db" -w /db/databases/bowtie2 "${BOWTIE2_IMG}" \
    bowtie2-build --large-index --ftabchars 4 --seed 42 --threads 1 -q /db/genomes.fna WoLr2

ls -l "${DB}/databases/bowtie2"

# Committed alignment of the fixture reads, made exactly as WOLTKA_ALIGN does
# (SHOGUN flags from config/modules.config, --no-head --no-unal, trimmed SAM)
# so WOLTKA_CLASSIFY / bin-script tests need no aligner. gzip -n keeps the
# file byte-stable across regenerations.
FIX="${REPO}/tests/data/woltka"
docker run --rm -u "$(id -u):$(id -g)" -v "${FIX}:/f" -w /f "${BOWTIE2_IMG}" bash -c '
    set -euo pipefail
    bowtie2 -p 1 -x db/databases/bowtie2/WoLr2 \
        -1 reads/wtest_R1.fastq.gz -2 reads/wtest_R2.fastq.gz -U reads/wtest_S.fastq.gz \
        --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" --seed 42 \
        --no-head --no-unal 2> /dev/null \
        | cut -f1-9 | sed "s/\$/\t*\t*/" | gzip -n > wtest.wol.sam.gz'
zcat "${FIX}/wtest.wol.sam.gz" | wc -l

# Per-sample WOLTKA_CLASSIFY outputs for the AGGREGATE_READ_FUNCTIONS module
# test: classify (the module's fixed flags) + bin/woltka_function_profile.py
# in the pinned woltka image, plus a versions.yml in the module's shape
WOLTKA_IMG="quay.io/biocontainers/woltka:0.1.7--pyhdfd78af_0"
docker run --rm -u "$(id -u):$(id -g)" -v "${FIX}:/f" -v "${REPO}/bin:/b:ro" -w /tmp "${WOLTKA_IMG}" sh -c '
    set -e
    woltka classify --input /f/wtest.wol.sam.gz --coords /f/db/proteins/coords.txt.xz \
        --no-demux --digits 6 --unassigned --to-tsv --output /tmp/orf.tsv > /dev/null
    cd /f && python3 /b/woltka_function_profile.py --profile /tmp/orf.tsv --db /f/db \
        --cog-def /f/cog.def.tab --sample-id wtest --prefix wtest'
printf '"WOLTKA_CLASSIFY":\n    woltka: 0.1.7\n    wol_db: WoLr2 test fixture (tests/bin/make_woltka_fixtures.py)\n    cog_def: cog.def.tab\n' \
    > "${FIX}/woltka_versions.yml"
