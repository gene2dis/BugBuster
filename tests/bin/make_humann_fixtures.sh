#!/bin/bash
#
# Regenerate the committed HUMAnN read-branch fixture (design doc T8c):
# tests/data/humann/{reads,db,htest_2_genefamilies.tsv,htest_3_reactions.tsv,
# htest_4_pathabundance.tsv,htest.humann_functions.tsv,htest.humann_summary.tsv,
# humann_versions.yml,htest.ags.tsv}. Everything derives from HUMAnN 4.0.0a2's
# bundled demo data inside the pinned HUMAnN image (see
# make_humann_fixtures.py for the layout and the synthetic KO/EC maps).
#
# The fixture HUMAnN run uses --bypass-prescreen (no MetaPhlAn database in CI),
# both subject-coverage thresholds at 0 (a 2,874-read subsample covers almost no
# gene to HUMAnN's default 50 %, leaving no pathways - the fixture tests the
# layout and wiring, not biology), --threads 1 and xipe off (the default) so
# the outputs are byte-stable. The same extra arguments are the HUMANN ext.args
# of the module nf-test (tests/modules/humann.config).
#
# Usage: tests/bin/make_humann_fixtures.sh   (needs python3 and docker)

set -euo pipefail

FIXTURE_ARGS="--bypass-prescreen --nucleotide-subject-coverage-threshold 0 --translated-subject-coverage-threshold 0"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMG="community.wave.seqera.io/library/python_metaphlan_diamond_bowtie2_pruned:386bc0e5651a4c44"
FIX="${REPO}/tests/data/humann"
STAGING="$(mktemp -d)"
trap 'rm -rf "${STAGING}"' EXIT

# 1. demo data out of the image (package data of the pinned humann 4.0.0a2)
docker run --rm -u "$(id -u):$(id -g)" -v "${STAGING}:/st" "${IMG}" bash -c '
    set -e
    d=/opt/conda/lib/python3.11/site-packages/humann
    cp ${d}/tests/data/demo.fastq /st/
    cp -r ${d}/data/chocophlan_DEMO /st/chocophlan
    cp -r ${d}/data/uniref_DEMO /st/uniref
    mkdir -p /st/utility
    for f in metacyc_reactions_level4ec_only.uniref.bz2 metacyc_pathways_structured_filtered_v24_subreactions \
             map_ko_name.txt.gz map_level4ec_name.txt.gz; do
        cp ${d}/data/utility_DEMO/${f} /st/utility/
    done'

mkdir -p "${FIX}"
python3 "${REPO}/tests/bin/make_humann_fixtures.py" stage1 "${STAGING}" "${FIX}"

# 2. a real HUMAnN run on the fixture, exactly as the HUMANN module runs it
#    (database flags, RPKs, temp removal) plus --bypass-prescreen
run_humann() {
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e FIXTURE_ARGS="${FIXTURE_ARGS}" -v "${FIX}:/f" -v "${REPO}/bin:/b:ro" -w /tmp "${IMG}" bash -c '
        set -e
        zcat /f/reads/htest_R1.fastq.gz /f/reads/htest_R2.fastq.gz > reads.fastq
        humann --input reads.fastq --input-format fastq --output out --output-basename htest \
            --threads 1 \
            --nucleotide-database /f/db/chocophlan --protein-database /f/db/uniref \
            --utility-database /f/db/utility_mapping \
            --pathways-database /f/db/utility_mapping/metacyc_reactions_level4ec_only.uniref.bz2,/f/db/utility_mapping/metacyc_pathways_structured_filtered_v24_subreactions \
            --count-normalization RPKs --remove-temp-output ${FIXTURE_ARGS} > /dev/null
        cp out/htest_2_genefamilies.tsv out/htest_3_reactions.tsv out/htest_4_pathabundance.tsv /f/
        echo $(( $(wc -l < reads.fastq) / 4 )) > /f/.input_reads'
}
run_humann

# 3. synthetic KO/EC maps over the families the run produced
python3 "${REPO}/tests/bin/make_humann_fixtures.py" stage2 "${FIX}" "${STAGING}"

# 4. committed composer outputs (inputs of the aggregation tests)
docker run --rm -u "$(id -u):$(id -g)" -v "${FIX}:/f" -v "${REPO}/bin:/b:ro" -w /f "${IMG}" \
    python3 /b/humann_function_profile.py \
        --genefamilies htest_2_genefamilies.tsv --pathabundance htest_4_pathabundance.tsv \
        --utility-db db/utility_mapping --basename htest --input-reads "$(cat "${FIX}/.input_reads")" \
        --sample-id htest --prefix htest
rm -f "${FIX}/.input_reads"

printf '"HUMANN":\n    humann: 4.0.0.alpha.2\n    metaphlan: 4.1.2\n    diamond: 2.1.24\n    bowtie2: 2.5.5\n    humann_db: fixture\n' \
    > "${FIX}/humann_versions.yml"
printf 'sample_id\taverage_genome_size_bp\tgenome_equivalents\ttotal_bases\nhtest\t3000000\t2.0\t6000000\n' \
    > "${FIX}/htest.ags.tsv"

ls -la "${FIX}"
