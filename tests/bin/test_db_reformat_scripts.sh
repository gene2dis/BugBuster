#!/bin/bash
#
# Regression tests for audit items #3 and #10 (plus the reformat scripts
# added since, e.g. sylph_db_reformat.sh):
#   - the pinned download container must contain every tool the bin/ DB
#     reformat scripts need (item #3), and
#   - the reformat scripts must fail loudly on partial/corrupt downloads
#     instead of caching a broken database with exit 0 (item #10).
#
# Each reformat script is run inside the same pinned image the format_db
# processes use, against tiny fixtures served over a local HTTP server.
#
# Requirements: docker, python3 (both present on GitHub ubuntu-latest runners).
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"

# Use the exact images pinned in the modules so the test cannot drift from them
DOWNLOAD_IMG=$(grep -o "'[^']*community.wave.seqera.io/library/wget:[^']*'" "${REPO_DIR}/modules/local/format_db/main.nf" \
    | tr -d "'" | grep -v '^oras://' | head -1)
# FORMAT_BAKTA_DB runs in the pinned bakta container (it needs
# amrfinder_update; the real binary is shadowed by a mock in these tests so
# no NCBI download happens - see MOCK_BIN below)
BAKTA_IMG=$(grep -o "'quay.io/biocontainers/bakta:[^']*'" "${REPO_DIR}/modules/local/format_db/main.nf" \
    | tr -d "'" | head -1)
# FORMAT_WOLTKA_DB runs in the pinned woltka container (it needs xz for the
# md5-of-uncompressed checks; its wget is busybox)
WOLTKA_IMG=$(grep -o "'quay.io/biocontainers/woltka:[^']*'" "${REPO_DIR}/modules/local/format_db/main.nf" \
    | tr -d "'" | head -1)
REPORT_IMG=$(grep -o "'[^']*community.wave.seqera.io/library/python_pandas[^']*'" "${REPO_DIR}/modules/local/taxonomy_report/main.nf" \
    | tr -d "'" | grep -v '^oras://' | head -1)
# FORMAT_SUPERFOCUS_DB runs in the pinned SUPER-FOCUS image (GNU wget + unzip)
SF_IMG=$(grep -o "'community.wave.seqera.io/library/super-focus[^']*'" "${REPO_DIR}/modules/local/format_db/main.nf" \
    | tr -d "'" | head -1)
# FORMAT_KMARVEL_DB runs in the pipeline-built K-MARVEL image (GNU wget +
# bzip2 for the CARD .tar.bz2)
KMARVEL_IMG=$(grep -o "'ghcr.io/gene2dis/bugbuster-kmarvel[^']*'" "${REPO_DIR}/modules/local/format_db/main.nf" \
    | tr -d "'" | head -1)

if [ -z "${DOWNLOAD_IMG}" ] || [ -z "${BAKTA_IMG}" ] || [ -z "${WOLTKA_IMG}" ] || [ -z "${REPORT_IMG}" ] || [ -z "${SF_IMG}" ] || [ -z "${KMARVEL_IMG}" ]; then
    echo "ERROR: could not extract pinned container images from the modules"
    exit 1
fi

echo "=== DB reformat script test suite ==="
echo "Download image: ${DOWNLOAD_IMG}"
echo "Bakta image:    ${BAKTA_IMG}"
echo "Woltka image:   ${WOLTKA_IMG}"
echo "Report image:   ${REPORT_IMG}"
echo "SUPER-FOCUS image: ${SF_IMG}"
echo "K-MARVEL image: ${KMARVEL_IMG}"
echo ""

TMP_DIR=$(mktemp -d)
SERVER_PID=""
cleanup() {
    [ -n "${SERVER_PID}" ] && kill "${SERVER_PID}" 2>/dev/null || true
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

FIXTURES="${TMP_DIR}/fixtures"
WORK="${TMP_DIR}/work"
mkdir -p "${FIXTURES}" "${WORK}"

#
# Build fixtures
#
build_fixtures() {
    local f="${FIXTURES}"
    local staging="${TMP_DIR}/staging"

    # Kraken2 DB tarball
    mkdir -p "${staging}/kraken"
    ( cd "${staging}/kraken" \
        && echo k > hash.k2d && echo k > opts.k2d && echo k > taxo.k2d \
        && tar -czf "${f}/k2_mini.tar.gz" hash.k2d opts.k2d taxo.k2d )

    # NCBI taxdump tarball
    mkdir -p "${staging}/taxdump"
    ( cd "${staging}/taxdump" \
        && echo n > nodes.dmp && echo n > names.dmp && echo n > division.dmp \
        && echo g > gc.prt && echo r > readme.txt \
        && tar -czf "${f}/taxdump.tar.gz" ./*.dmp gc.prt readme.txt )

    # CheckM2 DB tarball
    mkdir -p "${staging}/checkm2/CheckM2_database"
    ( cd "${staging}/checkm2" \
        && echo d > CheckM2_database/uniref100.KO.1.dmnd && echo c > CONTENTS.json \
        && tar -czf "${f}/checkm2_database.tar.gz" CheckM2_database CONTENTS.json )

    # GTDB-Tk data tarball (and an empty variant for the failure case)
    mkdir -p "${staging}/gtdbtk/release232/markers"
    ( cd "${staging}/gtdbtk" \
        && echo m > release232/markers/marker.txt \
        && tar -czf "${f}/gtdbtk_r232_data.tar.gz" release232 \
        && tar -czf "${f}/gtdbtk_dotslash.tar.gz" ./release232 )
    mkdir -p "${staging}/gtdbtk_empty/release232"
    ( cd "${staging}/gtdbtk_empty" && tar -czf "${f}/gtdbtk_empty.tar.gz" release232 )

    # Metabuli database tarball (refseq_standard layout: one directory with
    # the index files and the taxdump in ncbi-taxdump/), a variant without
    # the taxdump and one without db.parameters
    mkdir -p "${staging}/metabuli/refseq_standard/ncbi-taxdump"
    ( cd "${staging}/metabuli" \
        && for mf in db.parameters diffIdx info taxonomyDB; do echo m > "refseq_standard/${mf}"; done \
        && echo n > refseq_standard/ncbi-taxdump/nodes.dmp && echo n > refseq_standard/ncbi-taxdump/names.dmp \
        && tar -czf "${f}/metabuli_mini.tar.gz" refseq_standard \
        && rm -r refseq_standard/ncbi-taxdump \
        && tar -czf "${f}/metabuli_notaxdump.tar.gz" refseq_standard \
        && rm refseq_standard/db.parameters \
        && tar -czf "${f}/metabuli_noparams.tar.gz" refseq_standard )

    # geNomad database release tarball (one genomad_db/ directory) and one
    # whose version.txt names another release
    mkdir -p "${staging}/genomad/genomad_db"
    ( cd "${staging}/genomad" \
        && for gf in genomad_db genomad_db.dbtype genomad_marker_metadata.tsv; do echo g > "genomad_db/${gf}"; done \
        && echo 1.9 > genomad_db/version.txt \
        && tar -czf "${f}/genomad_db_v1.9.tar.gz" genomad_db \
        && echo 1.8 > genomad_db/version.txt \
        && tar -czf "${f}/genomad_db_v1.8.tar.gz" genomad_db )

    # Mock eggNOG emapper-3.0 data repository: the seven uncompressed data
    # files the reformat script fetches from the base URL
    mkdir -p "${f}/eggnogdb"
    for egf in eggnog.db eggnog.db.fieldpresence.bin eggnog.db.taxids.bin \
               eggnog.taxa.db eggnog.taxa.db.traverse.pkl eggnog_proteins.dmnd \
               go-basic.obo; do
        echo e > "${f}/eggnogdb/${egf}"
    done

    # Same eggNOG repo but with an empty eggnog.db (files are uncompressed, so
    # zero length is the detectable download corruption)
    cp -r "${f}/eggnogdb" "${f}/eggnogdb_empty"
    : > "${f}/eggnogdb_empty/eggnog.db"

    # Mock dbCAN release repository: the four CAZyme-annotation files the
    # reformat script fetches (note the remote dbCAN_sub.hmm underscore name)
    mkdir -p "${f}/dbcandb"
    for dbf in CAZy.dmnd dbCAN.hmm dbCAN_sub.hmm fam-substrate-mapping.tsv; do
        echo d > "${f}/dbcandb/${dbf}"
    done

    # Same dbCAN repo but with an empty CAZy.dmnd (uncompressed files, so
    # zero length is the detectable download corruption)
    cp -r "${f}/dbcandb" "${f}/dbcandb_empty"
    : > "${f}/dbcandb_empty/CAZy.dmnd"

    # Mock Bakta DB tarballs (.tar.xz, extracted with --strip-components=1):
    # full (db/) and light (db-light/) flavors, a wrong-schema variant, and
    # one whose bundled AMRFinderPlus DB is missing (truncated download)
    mkdir -p "${staging}/bakta_full/db/amrfinderplus-db"
    ( cd "${staging}/bakta_full" \
        && printf '{"date": "2025-02-24", "major": 6, "minor": 0, "type": "full"}\n' > db/version.json \
        && echo a > db/amrfinderplus-db/AMR.LIB && echo b > db/bakta.db \
        && tar -cJf "${f}/bakta_db_full.tar.xz" db )
    mkdir -p "${staging}/bakta_light/db-light/amrfinderplus-db"
    ( cd "${staging}/bakta_light" \
        && printf '{"date": "2025-02-24", "major": 6, "minor": 0, "type": "light"}\n' > db-light/version.json \
        && echo a > db-light/amrfinderplus-db/AMR.LIB && echo b > db-light/bakta.db \
        && tar -cJf "${f}/bakta_db_light.tar.xz" db-light )
    mkdir -p "${staging}/bakta_v5/db/amrfinderplus-db"
    ( cd "${staging}/bakta_v5" \
        && printf '{"date": "2023-02-20", "major": 5, "minor": 1, "type": "full"}\n' > db/version.json \
        && echo a > db/amrfinderplus-db/AMR.LIB && echo b > db/bakta.db \
        && tar -cJf "${f}/bakta_db_v5.tar.xz" db )
    mkdir -p "${staging}/bakta_noamr/db"
    ( cd "${staging}/bakta_noamr" \
        && printf '{"date": "2025-02-24", "major": 6, "minor": 0, "type": "full"}\n' > db/version.json \
        && echo b > db/bakta.db \
        && tar -cJf "${f}/bakta_db_noamr.tar.xz" db )

    # A corrupt (truncated) gzip tarball, and an xz twin for the Bakta script
    head -c 100 /dev/urandom > "${f}/corrupt.tar.gz"
    head -c 100 /dev/urandom > "${f}/corrupt.tar.xz"

    # sylph: the committed two-genome module fixture database and its
    # sylph-tax metadata, an empty database, a corrupt gzip metadata file and
    # a valid gzip that is not sylph-tax metadata
    cp "${REPO_DIR}/tests/data/sylph/sylph_test_db.syldb" "${f}/mini-c200.syldb"
    : > "${f}/empty-c200.syldb"
    cp "${REPO_DIR}/tests/data/sylph/sylph_tax_fixture.tsv.gz" "${f}/mini_metadata.tsv.gz"
    head -c 100 /dev/urandom > "${f}/corrupt_metadata.tsv.gz"
    printf 'just some text\n' | gzip -n > "${f}/plain_metadata.tsv.gz"
    # an HTML error page saved under each database name, a two-stage
    # database stand-in with the SY2D magic, and binary data without it
    printf '<!DOCTYPE html><html><body>503 Service Unavailable</body></html>\n' > "${f}/html-c200.syldb"
    cp "${f}/html-c200.syldb" "${f}/html-c200.syl2db"
    { printf 'SY2D'; head -c 2048 /dev/zero; } > "${f}/mini-c200.syl2db"
    { printf 'XXXX'; head -c 2048 /dev/zero; } > "${f}/nomagic-c200.syl2db"
}
build_fixtures

# Woltka WoLr2-layout fixtures, copied from the committed read-branch fixture
# (tests/data/woltka/db, made by tests/bin/make_woltka_fixtures.sh)
build_wol_fixture() {
    local dest="${FIXTURES}/$1"
    mkdir -p "${dest}"
    cp -r "${REPO_DIR}/tests/data/woltka/db/databases" "${REPO_DIR}/tests/data/woltka/db/proteins" \
        "${REPO_DIR}/tests/data/woltka/db/function" "${dest}/"
}
# NCBI COG definitions table fixtures (cog-24.def.tab layout + the release's
# checksums.md5 next to it), plus a bad-md5 and a wrong-layout variant
for v in cog cog_badmd5 cog_badlayout; do mkdir -p "${FIXTURES}/${v}"; done
cp "${REPO_DIR}/tests/data/functional/cog-fixture.def.tab" "${FIXTURES}/cog/cog-24.def.tab"
cp "${FIXTURES}/cog/cog-24.def.tab" "${FIXTURES}/cog_badmd5/cog-24.def.tab"
printf 'not a cog table\n' > "${FIXTURES}/cog_badlayout/cog-24.def.tab"
for v in cog cog_badlayout; do
    (cd "${FIXTURES}/${v}" && md5sum cog-24.def.tab > checksums.md5)
done
echo "00000000000000000000000000000000  cog-24.def.tab" > "${FIXTURES}/cog_badmd5/checksums.md5"

build_wol_fixture wol2
build_wol_fixture wol2_badmd5
echo "00000000000000000000000000000000" > "${FIXTURES}/wol2_badmd5/proteins/coords.txt.md5"
build_wol_fixture wol2_emptyindex
: > "${FIXTURES}/wol2_emptyindex/databases/bowtie2/WoLr2.3.bt2l"
build_wol_fixture wol2_missing
rm "${FIXTURES}/wol2_missing/function/metacyc/pathway_name.txt"

# SUPER-FOCUS DB_90 archives as figshare serves them, zipped from the
# committed read-branch fixture (tests/data/superfocus/sf_db, made by
# tests/bin/make_superfocus_fixtures.sh): diamond .dmnd at the archive root,
# mmseqs2 files inside a folder (the script unpacks with -j), plus broken
# variants, and the database_PKs.txt table
mkdir -p "${FIXTURES}/sf"
python3 - "${REPO_DIR}/tests/data/superfocus/sf_db/db" "${FIXTURES}/sf" <<'PY'
import pathlib, sys, zipfile
db, out = map(pathlib.Path, sys.argv[1:])
mmseqs = sorted((db / 'static' / 'mmseqs2').iterdir())
with zipfile.ZipFile(out / 'diamond90.zip', 'w') as z:
    z.write(db / 'static' / 'diamond' / '90_clusters.db.dmnd', '90_clusters.db.dmnd')
with zipfile.ZipFile(out / 'mmseqs90.zip', 'w') as z:
    for f in mmseqs:
        z.write(f, f'mmseqs_90/{f.name}')
with zipfile.ZipFile(out / 'mmseqs90_nodbtype.zip', 'w') as z:
    for f in mmseqs:
        if f.name != '90_clusters.db.dbtype':
            z.write(f, f'mmseqs_90/{f.name}')
with zipfile.ZipFile(out / 'wrong.zip', 'w') as z:
    z.writestr('README.txt', 'not a SUPER-FOCUS database\n')
PY
cp "${REPO_DIR}/tests/data/superfocus/sf_db/db/database_PKs.txt" "${FIXTURES}/sf/database_PKs.txt"
printf 'not\tthe\tsubsystem table\n' > "${FIXTURES}/sf/bad_PKs.txt"

# HUMAnN database archives as huttenhower.sph.harvard.edu / the MetaPhlAn
# server publish them, built from the committed read-branch fixture database
# (tests/data/humann/db, made by tests/bin/make_humann_fixtures.sh): flat
# tar.gz archives (one ChocoPhlAn variant wrapped in a folder), plain MetaPhlAn
# tars with their published-style .md5 files (the Bowtie2 index inside a
# folder), plus broken variants
HU="${FIXTURES}/hu"
HDB="${REPO_DIR}/tests/data/humann/db"
mkdir -p "${HU}/wrap/chocophlan_v4" "${HU}/mpa" "${HU}/bt2/bowtie2_indexes"
tar -czf "${HU}/chocophlan.tar.gz" -C "${HDB}/chocophlan" .
cp "${HDB}"/chocophlan/* "${HU}/wrap/chocophlan_v4/"
tar -czf "${HU}/chocophlan_wrapped.tar.gz" -C "${HU}/wrap" chocophlan_v4
tar -czf "${HU}/uniref.tar.gz" -C "${HDB}/uniref" .
tar -czf "${HU}/utility.tar.gz" -C "${HDB}/utility_mapping" .
tar -czf "${HU}/utility_noko.tar.gz" -C "${HDB}/utility_mapping" --exclude map_ko_uniref90.txt.gz .
tar -czf "${HU}/empty.tar.gz" -C "${HU}/wrap" --files-from /dev/null
head -c 2000 "${HU}/chocophlan.tar.gz" > "${HU}/truncated.tar.gz"
IDX=mpa_vOct22_CHOCOPhlAnSGB_202403
printf 'pkl fixture\n' > "${HU}/mpa/${IDX}.pkl"
tar -cf "${HU}/${IDX}.tar" -C "${HU}/mpa" "${IDX}.pkl"
for ext in 1 2 3 4 rev.1 rev.2; do printf 'bt2l %s\n' "${ext}" > "${HU}/bt2/bowtie2_indexes/${IDX}.${ext}.bt2l"; done
tar -cf "${HU}/${IDX}_bt2.tar" -C "${HU}/bt2" bowtie2_indexes
rm "${HU}/bt2/bowtie2_indexes/${IDX}.rev.2.bt2l"
tar -cf "${HU}/bt2_incomplete.tar" -C "${HU}/bt2" bowtie2_indexes
( cd "${HU}" && md5sum "${IDX}.tar" > "${IDX}.md5" && md5sum "${IDX}_bt2.tar" > "${IDX}_bt2.md5" )
printf '%s  %s\n' "$(md5sum "${HU}/bt2_incomplete.tar" | cut -d' ' -f1)" "${IDX}_bt2.tar" > "${HU}/bt2_incomplete.md5"
printf '00000000000000000000000000000000  %s\n' "${IDX}_bt2.tar" > "${HU}/bad_bt2.md5"

# Mock CARD data tarballs (.tar.bz2, files at the archive root with a
# leading ./ as CARD ships them): complete, and one without the protein
# variant model FASTA
KM="${FIXTURES}/km"
mkdir -p "${KM}/card" "${KM}/card_novar"
printf '>gb|A|ARO:1|blaX\nMKLV\n' > "${KM}/card/protein_fasta_protein_homolog_model.fasta"
printf '>gb|B|ARO:2|ompX\nMQQA\n' > "${KM}/card/protein_fasta_protein_overexpression_model.fasta"
printf '>gb|C|ARO:3|gyrA\nMSDL\n' > "${KM}/card/protein_fasta_protein_variant_model.fasta"
printf 'ARO Accession\tCVTERM ID\nARO:1\t1\n' > "${KM}/card/aro_index.tsv"
printf 'Accession\tName\tModel Type\n3\tgyrA\tprotein variant model\n' > "${KM}/card/snps.txt"
printf 'nucleotide\n' > "${KM}/card/nucleotide_fasta_protein_homolog_model.fasta"
tar -cjf "${KM}/card.tar.bz2" -C "${KM}/card" .
cp "${KM}"/card/* "${KM}/card_novar/"
rm "${KM}/card_novar/protein_fasta_protein_variant_model.fasta"
tar -cjf "${KM}/card_novar.tar.bz2" -C "${KM}/card_novar" .

#
# Local HTTP server for the fixtures
#
PORT=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1]); s.close()')
python3 -m http.server --bind 127.0.0.1 --directory "${FIXTURES}" "${PORT}" > /dev/null 2>&1 &
SERVER_PID=$!
for _ in $(seq 1 50); do
    curl -sf "http://127.0.0.1:${PORT}/" > /dev/null && break
    sleep 0.1
done
BASE_URL="http://127.0.0.1:${PORT}"

#
# Helpers
#
PASS=0
FAIL=0

run_script() {
    # run_script <workdir> <script-name> [args...] — inside the download image
    # (or RUN_IMG when set, for scripts pinned to a different container).
    # When MOCK_BIN is set, each executable in it is bind-mounted OVER the
    # container's real /usr/local/bin/<name> (used to stub out
    # amrfinder_update's NCBI download in the bakta tests; a PATH override
    # does not survive the biocontainer entrypoint's PATH rewriting).
    local workdir="$1" script="$2"
    shift 2
    local mock_args=()
    if [ -n "${MOCK_BIN:-}" ]; then
        local mock_file
        for mock_file in "${MOCK_BIN}"/*; do
            mock_args+=(-v "${mock_file}:/usr/local/bin/$(basename "${mock_file}"):ro")
        done
    fi
    docker run --rm --network host -u "$(id -u):$(id -g)" -e HOME=/tmp \
        -e WOLTKA_DB_RETRIES=1 -e WOLTKA_DB_RETRY_WAIT=0 \
        -e HUMANN_DB_RETRIES=1 -e HUMANN_DB_RETRY_WAIT=0 \
        "${mock_args[@]}" \
        -v "${REPO_DIR}/bin:/pipeline_bin:ro" \
        -v "${workdir}:/dbwork" -w /dbwork \
        "${RUN_IMG:-${DOWNLOAD_IMG}}" bash "/pipeline_bin/${script}" "$@"
}

expect_pass() {
    local desc="$1"; shift
    local workdir="${WORK}/$(echo "${desc}" | tr ' /:' '___')"
    mkdir -p "${workdir}"
    if run_script "${workdir}" "$@" > "${workdir}.log" 2>&1; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
        LAST_WORKDIR="${workdir}"
    else
        echo "✗ ${desc} — expected success, got failure:"
        tail -5 "${workdir}.log" | sed 's/^/    /'
        FAIL=$((FAIL + 1))
        LAST_WORKDIR="${workdir}"
    fi
}

expect_fail() {
    local desc="$1"; shift
    local workdir="${WORK}/$(echo "${desc}" | tr ' /:' '___')"
    mkdir -p "${workdir}"
    if run_script "${workdir}" "$@" > "${workdir}.log" 2>&1; then
        echo "✗ ${desc} — expected nonzero exit, script succeeded"
        FAIL=$((FAIL + 1))
    else
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    fi
}

check_file() {
    # check_file <desc> <path-relative-to-last-workdir>
    if [ -e "${LAST_WORKDIR}/$2" ]; then
        echo "✓ $1"
        PASS=$((PASS + 1))
    else
        echo "✗ $1 — missing: $2"
        FAIL=$((FAIL + 1))
    fi
}

#
# kraken_db_reformat.sh
#
echo "--- kraken_db_reformat.sh ---"
expect_pass "kraken: valid tarball" kraken_db_reformat.sh "${BASE_URL}/k2_mini.tar.gz"
check_file "kraken: hash.k2d extracted" "k2_ref_db/hash.k2d"
check_file "kraken: taxo.k2d extracted" "k2_ref_db/taxo.k2d"
expect_fail "kraken: corrupt tarball fails" kraken_db_reformat.sh "${BASE_URL}/corrupt.tar.gz"
expect_fail "kraken: 404 URL fails" kraken_db_reformat.sh "${BASE_URL}/no_such_file.tar.gz"
expect_fail "kraken: tarball without k2d files fails" kraken_db_reformat.sh "${BASE_URL}/taxdump.tar.gz"

#
# sylph_db_reformat.sh
#
echo "--- sylph_db_reformat.sh ---"
expect_pass "sylph: valid database" sylph_db_reformat.sh db "${BASE_URL}/mini-c200.syldb"
check_file "sylph: database saved under its URL name" "sylph_db/mini-c200.syldb"
expect_pass "sylph: valid taxonomy metadata" sylph_db_reformat.sh taxonomy "${BASE_URL}/mini_metadata.tsv.gz"
check_file "sylph: metadata saved under its URL name" "sylph_db/mini_metadata.tsv.gz"
expect_fail "sylph: empty database fails" sylph_db_reformat.sh db "${BASE_URL}/empty-c200.syldb"
expect_pass "sylph: two-stage database with SY2D magic" sylph_db_reformat.sh db "${BASE_URL}/mini-c200.syl2db"
expect_fail "sylph: HTML page saved as .syldb fails" sylph_db_reformat.sh db "${BASE_URL}/html-c200.syldb"
expect_fail "sylph: HTML page saved as .syl2db fails" sylph_db_reformat.sh db "${BASE_URL}/html-c200.syl2db"
expect_fail "sylph: .syl2db without SY2D magic fails" sylph_db_reformat.sh db "${BASE_URL}/nomagic-c200.syl2db"
expect_fail "sylph: 404 database URL fails" sylph_db_reformat.sh db "${BASE_URL}/no_such.syl2db"
expect_fail "sylph: non-sylph database name fails" sylph_db_reformat.sh db "${BASE_URL}/k2_mini.tar.gz"
expect_fail "sylph: corrupt metadata gzip fails" sylph_db_reformat.sh taxonomy "${BASE_URL}/corrupt_metadata.tsv.gz"
expect_fail "sylph: gzip without accession/lineage layout fails" sylph_db_reformat.sh taxonomy "${BASE_URL}/plain_metadata.tsv.gz"
expect_fail "sylph: unknown mode fails" sylph_db_reformat.sh index "${BASE_URL}/mini-c200.syldb"

#
# checkm2_db_reformat.sh
#
echo "--- checkm2_db_reformat.sh ---"
expect_pass "checkm2: valid tarball" checkm2_db_reformat.sh "${BASE_URL}/checkm2_database.tar.gz"
check_file "checkm2: dmnd file extracted" "uniref100.KO.1.dmnd"
expect_fail "checkm2: corrupt tarball fails" checkm2_db_reformat.sh "${BASE_URL}/corrupt.tar.gz"
expect_fail "checkm2: tarball without dmnd fails" checkm2_db_reformat.sh "${BASE_URL}/taxdump.tar.gz"

#
# gtdb-tk_db_reformat.sh
#
echo "--- gtdb-tk_db_reformat.sh ---"
expect_pass "gtdbtk: valid tarball" gtdb-tk_db_reformat.sh "${BASE_URL}/gtdbtk_r232_data.tar.gz"
check_file "gtdbtk: release files extracted" "gtdbtk_db/markers/marker.txt"
expect_fail "gtdbtk: corrupt tarball fails" gtdb-tk_db_reformat.sh "${BASE_URL}/corrupt.tar.gz"
expect_fail "gtdbtk: empty package fails" gtdb-tk_db_reformat.sh "${BASE_URL}/gtdbtk_empty.tar.gz"
expect_pass "gtdbtk: ./release232 entries" gtdb-tk_db_reformat.sh "${BASE_URL}/gtdbtk_dotslash.tar.gz"
check_file "gtdbtk: ./ prefix still gives gtdbtk_db/markers" "gtdbtk_db/markers/marker.txt"

#
# metabuli_db_reformat.sh
#
echo "--- metabuli_db_reformat.sh ---"
expect_pass "metabuli: valid tarball" metabuli_db_reformat.sh "${BASE_URL}/metabuli_mini.tar.gz"
check_file "metabuli: index files in metabuli_db" "metabuli_db/diffIdx"
check_file "metabuli: taxdump kept in metabuli_db" "metabuli_db/ncbi-taxdump/nodes.dmp"
expect_fail "metabuli: corrupt tarball fails" metabuli_db_reformat.sh "${BASE_URL}/corrupt.tar.gz"
expect_fail "metabuli: tarball without db.parameters fails" metabuli_db_reformat.sh "${BASE_URL}/metabuli_noparams.tar.gz"
expect_fail "metabuli: database without a taxdump fails" metabuli_db_reformat.sh "${BASE_URL}/metabuli_notaxdump.tar.gz"
expect_fail "metabuli: 404 URL fails" metabuli_db_reformat.sh "${BASE_URL}/no_such_file.tar.gz"

#
# genomad_db_reformat.sh
#
echo "--- genomad_db_reformat.sh ---"
GENOMAD_MD5=$(md5sum "${FIXTURES}/genomad_db_v1.9.tar.gz" | cut -d' ' -f1)
GENOMAD_OLD_MD5=$(md5sum "${FIXTURES}/genomad_db_v1.8.tar.gz" | cut -d' ' -f1)
expect_pass "genomad: valid release" genomad_db_reformat.sh "${BASE_URL}/genomad_db_v1.9.tar.gz" "${GENOMAD_MD5}" 1.9
check_file "genomad: database files in genomad_db" "genomad_db/genomad_db"
check_file "genomad: version.txt kept" "genomad_db/version.txt"
expect_fail "genomad: wrong md5 fails" genomad_db_reformat.sh "${BASE_URL}/genomad_db_v1.9.tar.gz" \
    00000000000000000000000000000000 1.9
expect_fail "genomad: another release than the registry's fails" genomad_db_reformat.sh \
    "${BASE_URL}/genomad_db_v1.8.tar.gz" "${GENOMAD_OLD_MD5}" 1.9
expect_fail "genomad: 404 URL fails" genomad_db_reformat.sh "${BASE_URL}/no_such_file.tar.gz" "${GENOMAD_MD5}" 1.9

#
# eggnog_db_reformat.sh
#
echo "--- eggnog_db_reformat.sh ---"
expect_pass "eggnog: valid data files" eggnog_db_reformat.sh "${BASE_URL}/eggnogdb"
check_file "eggnog: eggnog.db in eggnog_db" "eggnog_db/eggnog.db"
check_file "eggnog: taxid cache in eggnog_db" "eggnog_db/eggnog.db.taxids.bin"
check_file "eggnog: eggnog.taxa.db in eggnog_db" "eggnog_db/eggnog.taxa.db"
check_file "eggnog: diamond db in eggnog_db" "eggnog_db/eggnog_proteins.dmnd"
check_file "eggnog: go-basic.obo in eggnog_db" "eggnog_db/go-basic.obo"
expect_fail "eggnog: empty eggnog.db fails" eggnog_db_reformat.sh "${BASE_URL}/eggnogdb_empty"
expect_fail "eggnog: missing files fail" eggnog_db_reformat.sh "${BASE_URL}/no_such_dir"

#
# dbcan_db_reformat.sh
#
echo "--- dbcan_db_reformat.sh ---"
expect_pass "dbcan: valid data files" dbcan_db_reformat.sh "${BASE_URL}/dbcandb"
check_file "dbcan: CAZy.dmnd in dbcan_db" "dbcan_db/CAZy.dmnd"
check_file "dbcan: dbCAN.hmm in dbcan_db" "dbcan_db/dbCAN.hmm"
check_file "dbcan: dbCAN_sub.hmm saved as dbCAN-sub.hmm" "dbcan_db/dbCAN-sub.hmm"
check_file "dbcan: substrate mapping in dbcan_db" "dbcan_db/fam-substrate-mapping.tsv"
check_file "dbcan: DB_VERSION provenance file written" "dbcan_db/DB_VERSION"
expect_fail "dbcan: empty CAZy.dmnd fails" dbcan_db_reformat.sh "${BASE_URL}/dbcandb_empty"
expect_fail "dbcan: missing files fail" dbcan_db_reformat.sh "${BASE_URL}/no_such_dir"

#
# bakta_db_reformat.sh (runs in the bakta image FORMAT_BAKTA_DB pins; the
# real amrfinder_update is shadowed by a mock so the test never hits NCBI,
# while the script's update invocation and post-update verification still run)
#
echo "--- bakta_db_reformat.sh ---"
RUN_IMG="${BAKTA_IMG}"
MOCK_OK="${TMP_DIR}/mock_ok"
MOCK_FAIL="${TMP_DIR}/mock_fail"
mkdir -p "${MOCK_OK}" "${MOCK_FAIL}"
cat > "${MOCK_OK}/amrfinder_update" <<'EOF'
#!/bin/bash
# Mock: emulate a successful AMRFinderPlus DB refresh without the NCBI download
db=""
while [ $# -gt 0 ]; do
    case "$1" in
        --database) db="$2"; shift 2 ;;
        *) shift ;;
    esac
done
[ -n "$db" ] || { echo "mock amrfinder_update: no --database" >&2; exit 2; }
mkdir -p "$db/2099-01-01.1"
echo "2099-01-01.1" > "$db/2099-01-01.1/version.txt"
ln -sfn 2099-01-01.1 "$db/latest"
EOF
printf '#!/bin/bash\necho "mock amrfinder_update: simulated download failure" >&2\nexit 1\n' > "${MOCK_FAIL}/amrfinder_update"
chmod +x "${MOCK_OK}/amrfinder_update" "${MOCK_FAIL}/amrfinder_update"

MOCK_BIN="${MOCK_OK}"
expect_pass "bakta: valid full tarball" bakta_db_reformat.sh "${BASE_URL}/bakta_db_full.tar.xz"
check_file "bakta: version.json extracted (top dir stripped)" "bakta_db/version.json"
check_file "bakta: bundled amrfinderplus-db present" "bakta_db/amrfinderplus-db/AMR.LIB"
check_file "bakta: DB_VERSION provenance file written" "bakta_db/DB_VERSION"
if grep -q '^db-full ' "${LAST_WORKDIR}/bakta_db/DB_VERSION" 2>/dev/null; then
    echo "✓ bakta: DB_VERSION records the full flavor"
    PASS=$((PASS + 1))
else
    echo "✗ bakta: DB_VERSION does not record the full flavor"
    FAIL=$((FAIL + 1))
fi
if grep -q 'amrfinderplus-db 2099-01-01.1' "${LAST_WORKDIR}/bakta_db/DB_VERSION" 2>/dev/null \
    && [ "$(wc -l < "${LAST_WORKDIR}/bakta_db/DB_VERSION")" -eq 1 ]; then
    echo "✓ bakta: DB_VERSION records the refreshed AMRFinderPlus version on one line"
    PASS=$((PASS + 1))
else
    echo "✗ bakta: DB_VERSION missing the AMRFinderPlus version or not single-line"
    FAIL=$((FAIL + 1))
fi
check_file "bakta: refreshed amrfinderplus-db latest present" "bakta_db/amrfinderplus-db/latest/version.txt"
expect_pass "bakta: valid light tarball (db-light/ top dir)" bakta_db_reformat.sh "${BASE_URL}/bakta_db_light.tar.xz"
if grep -q '^db-light ' "${LAST_WORKDIR}/bakta_db/DB_VERSION" 2>/dev/null; then
    echo "✓ bakta: DB_VERSION records the light flavor"
    PASS=$((PASS + 1))
else
    echo "✗ bakta: DB_VERSION does not record the light flavor"
    FAIL=$((FAIL + 1))
fi
expect_fail "bakta: schema-5 version.json fails" bakta_db_reformat.sh "${BASE_URL}/bakta_db_v5.tar.xz"
expect_fail "bakta: missing amrfinderplus-db fails" bakta_db_reformat.sh "${BASE_URL}/bakta_db_noamr.tar.xz"
expect_fail "bakta: corrupt tarball fails" bakta_db_reformat.sh "${BASE_URL}/corrupt.tar.xz"
expect_fail "bakta: 404 URL fails" bakta_db_reformat.sh "${BASE_URL}/no_such_file.tar.xz"
MOCK_BIN="${MOCK_FAIL}"
expect_fail "bakta: amrfinder_update failure fails the provisioning" bakta_db_reformat.sh "${BASE_URL}/bakta_db_full.tar.xz"
MOCK_BIN=""
RUN_IMG=""

#
# cog_db_reformat.sh (shared download image, FORMAT_COG_DB)
#
echo "--- cog_db_reformat.sh ---"
expect_pass "cog: valid table + matching checksums.md5" cog_db_reformat.sh "${BASE_URL}/cog/cog-24.def.tab"
check_file "cog: table kept under its release name" "cog-24.def.tab"
if [ ! -e "${LAST_WORKDIR}/checksums.md5" ]; then
    echo "✓ cog: checksums.md5 removed after verification"; PASS=$((PASS + 1))
else
    echo "✗ cog: checksums.md5 left in the output dir"; FAIL=$((FAIL + 1))
fi
expect_fail "cog: md5 mismatch fails" cog_db_reformat.sh "${BASE_URL}/cog_badmd5/cog-24.def.tab"
expect_fail "cog: non-COG layout fails" cog_db_reformat.sh "${BASE_URL}/cog_badlayout/cog-24.def.tab"
expect_fail "cog: 404 URL fails" cog_db_reformat.sh "${BASE_URL}/no_such_dir/cog-24.def.tab"

#
# woltka_db_reformat.sh (runs in the woltka image FORMAT_WOLTKA_DB pins)
#
echo "--- woltka_db_reformat.sh ---"
RUN_IMG="${WOLTKA_IMG}"
expect_pass "woltka: valid WoLr2-layout files" woltka_db_reformat.sh "${BASE_URL}/wol2"
check_file "woltka: Bowtie2 index in the WoLr2 layout" "woltka_db/databases/bowtie2/WoLr2.rev.1.bt2l"
check_file "woltka: ORF coordinates" "woltka_db/proteins/coords.txt.xz"
check_file "woltka: ORF lengths" "woltka_db/proteins/length.map.xz"
check_file "woltka: KEGG ORF map" "woltka_db/function/kegg/orf-to-ko.map.xz"
check_file "woltka: KEGG KO-to-EC map" "woltka_db/function/kegg/ko-to-ec.map"
check_file "woltka: MetaCyc reaction-to-pathway map" "woltka_db/function/metacyc/reaction-to-pathway.map"
check_file "woltka: Pfam ORF map" "woltka_db/function/pfam/orf-to-pfam.map.xz"
if grep -q '^WoLr2 (' "${LAST_WORKDIR}/woltka_db/DB_VERSION" 2>/dev/null; then
    echo "✓ woltka: DB_VERSION records the WoLr2 release"
    PASS=$((PASS + 1))
else
    echo "✗ woltka: DB_VERSION missing or does not record the release"
    FAIL=$((FAIL + 1))
fi
expect_fail "woltka: md5 mismatch (uncompressed content) fails" woltka_db_reformat.sh "${BASE_URL}/wol2_badmd5"
expect_fail "woltka: empty index file fails" woltka_db_reformat.sh "${BASE_URL}/wol2_emptyindex"
expect_fail "woltka: missing map file fails" woltka_db_reformat.sh "${BASE_URL}/wol2_missing"
expect_fail "woltka: 404 base URL fails" woltka_db_reformat.sh "${BASE_URL}/no_such_dir"
RUN_IMG=""

#
# superfocus_db_reformat.sh
#
echo "--- superfocus_db_reformat.sh ---"
RUN_IMG="${SF_IMG}"
sf_md5() { md5sum "${FIXTURES}/sf/$1" | cut -d' ' -f1; }
SF_PKS=("${BASE_URL}/sf/database_PKs.txt" "$(sf_md5 database_PKs.txt)")
expect_pass "superfocus: diamond DB_90 archive" superfocus_db_reformat.sh diamond \
    "${BASE_URL}/sf/diamond90.zip" "$(sf_md5 diamond90.zip)" "${SF_PKS[@]}" "SUPER-FOCUS test"
check_file "superfocus: database_PKs.txt in the root's db/" "superfocus_db/db/database_PKs.txt"
check_file "superfocus: diamond .dmnd in db/static/diamond/" "superfocus_db/db/static/diamond/90_clusters.db.dmnd"
if grep -q '^SUPER-FOCUS test; diamond DB_90 (' "${LAST_WORKDIR}/superfocus_db/DB_VERSION" 2>/dev/null \
        && [ "$(wc -l < "${LAST_WORKDIR}/superfocus_db/DB_VERSION")" -eq 1 ]; then
    echo "✓ superfocus: single-line DB_VERSION records release, aligner and archive"
    PASS=$((PASS + 1))
else
    echo "✗ superfocus: DB_VERSION missing, multi-line or incomplete"
    FAIL=$((FAIL + 1))
fi
if [ -z "$(ls "${LAST_WORKDIR}" | grep -v '^superfocus_db$')" ]; then
    echo "✓ superfocus: downloaded archive removed after unpacking"
    PASS=$((PASS + 1))
else
    echo "✗ superfocus: leftover files next to superfocus_db/"
    FAIL=$((FAIL + 1))
fi
expect_pass "superfocus: mmseqs2 DB_90 archive (folder inside the zip)" superfocus_db_reformat.sh mmseqs2 \
    "${BASE_URL}/sf/mmseqs90.zip" "$(sf_md5 mmseqs90.zip)" "${SF_PKS[@]}" "SUPER-FOCUS test"
check_file "superfocus: mmseqs2 database unpacked flat" "superfocus_db/db/static/mmseqs2/90_clusters.db"
check_file "superfocus: mmseqs2 .dbtype present" "superfocus_db/db/static/mmseqs2/90_clusters.db.dbtype"
expect_fail "superfocus: archive md5 mismatch fails" superfocus_db_reformat.sh diamond \
    "${BASE_URL}/sf/diamond90.zip" 00000000000000000000000000000000 "${SF_PKS[@]}" "SUPER-FOCUS test"
expect_fail "superfocus: database_PKs.txt md5 mismatch fails" superfocus_db_reformat.sh diamond \
    "${BASE_URL}/sf/diamond90.zip" "$(sf_md5 diamond90.zip)" "${BASE_URL}/sf/database_PKs.txt" 00000000000000000000000000000000 "SUPER-FOCUS test"
expect_fail "superfocus: database_PKs.txt without the subsystem header fails" superfocus_db_reformat.sh diamond \
    "${BASE_URL}/sf/diamond90.zip" "$(sf_md5 diamond90.zip)" "${BASE_URL}/sf/bad_PKs.txt" "$(sf_md5 bad_PKs.txt)" "SUPER-FOCUS test"
expect_fail "superfocus: archive without the aligner's database fails" superfocus_db_reformat.sh diamond \
    "${BASE_URL}/sf/wrong.zip" "$(sf_md5 wrong.zip)" "${SF_PKS[@]}" "SUPER-FOCUS test"
expect_fail "superfocus: diamond archive given for mmseqs2 fails" superfocus_db_reformat.sh mmseqs2 \
    "${BASE_URL}/sf/diamond90.zip" "$(sf_md5 diamond90.zip)" "${SF_PKS[@]}" "SUPER-FOCUS test"
expect_fail "superfocus: incomplete mmseqs2 database (no .dbtype) fails" superfocus_db_reformat.sh mmseqs2 \
    "${BASE_URL}/sf/mmseqs90_nodbtype.zip" "$(sf_md5 mmseqs90_nodbtype.zip)" "${SF_PKS[@]}" "SUPER-FOCUS test"
expect_fail "superfocus: unknown aligner fails" superfocus_db_reformat.sh rapsearch \
    "${BASE_URL}/sf/diamond90.zip" "$(sf_md5 diamond90.zip)" "${SF_PKS[@]}" "SUPER-FOCUS test"
expect_fail "superfocus: 404 archive URL fails" superfocus_db_reformat.sh diamond \
    "${BASE_URL}/sf/no_such.zip" 00000000000000000000000000000000 "${SF_PKS[@]}" "SUPER-FOCUS test"
RUN_IMG=""

#
# kmarvel_db_reformat.sh
#
echo "--- kmarvel_db_reformat.sh ---"
RUN_IMG="${KMARVEL_IMG}"
km_md5() { md5sum "${FIXTURES}/km/$1" | cut -d' ' -f1; }
expect_pass "kmarvel: CARD data tarball" kmarvel_db_reformat.sh \
    "${BASE_URL}/km/card.tar.bz2" "$(km_md5 card.tar.bz2)" "CARD test"
check_file "kmarvel: metadata.tsv (aro_index.tsv) in kmarvel_db" "kmarvel_db/metadata.tsv"
check_file "kmarvel: snps.txt in kmarvel_db" "kmarvel_db/snps.txt"
if [ "$(grep -c '^>' "${LAST_WORKDIR}/kmarvel_db/CARD_protein_all.fasta" 2>/dev/null)" = "3" ] \
        && [ "$(cat "${LAST_WORKDIR}/kmarvel_db/DB_VERSION" 2>/dev/null)" = "CARD test" ]; then
    echo "✓ kmarvel: homolog + overexpression + variant FASTAs concatenated, DB_VERSION written"
    PASS=$((PASS + 1))
else
    echo "✗ kmarvel: CARD_protein_all.fasta lacks the 3 model FASTAs, or DB_VERSION is wrong"
    FAIL=$((FAIL + 1))
fi
if [ -z "$(ls "${LAST_WORKDIR}" | grep -v '^kmarvel_db$')" ] \
        && [ "$(LC_ALL=C ls "${LAST_WORKDIR}/kmarvel_db" | tr '\n' ' ')" = "CARD_protein_all.fasta CARD_protein_all.fasta.pdb CARD_protein_all.fasta.phr CARD_protein_all.fasta.pin CARD_protein_all.fasta.pjs CARD_protein_all.fasta.pot CARD_protein_all.fasta.psq CARD_protein_all.fasta.ptf CARD_protein_all.fasta.pto DB_VERSION metadata.tsv snps.txt " ]; then
    echo "✓ kmarvel: only the four database files and the FASTA's BLAST index are left"
    PASS=$((PASS + 1))
else
    echo "✗ kmarvel: leftover files next to or inside kmarvel_db/"
    FAIL=$((FAIL + 1))
fi
expect_fail "kmarvel: tarball md5 mismatch fails" kmarvel_db_reformat.sh \
    "${BASE_URL}/km/card.tar.bz2" 00000000000000000000000000000000 "CARD test"
expect_fail "kmarvel: tarball without the variant model FASTA fails" kmarvel_db_reformat.sh \
    "${BASE_URL}/km/card_novar.tar.bz2" "$(km_md5 card_novar.tar.bz2)" "CARD test"
expect_fail "kmarvel: 404 tarball URL fails" kmarvel_db_reformat.sh \
    "${BASE_URL}/km/no_such.tar.bz2" 00000000000000000000000000000000 "CARD test"
RUN_IMG=""

#
# humann_db_reformat.sh
#
echo "--- humann_db_reformat.sh ---"
H="${BASE_URL}/hu"
IDX=mpa_vOct22_CHOCOPhlAnSGB_202403
# hu_args <chocophlan> <uniref> <utility> <bt2-tar-url> <bt2-md5-url>
hu_args() {
    echo v4_alpha-test "${H}/$1" "${H}/$2" "${H}/$3" "${IDX}" "${H}/${IDX}.tar" "${H}/${IDX}.md5" "${H}/$4" "${H}/$5" "HUMAnN-test"
}
expect_pass "humann: full database set" humann_db_reformat.sh \
    $(hu_args chocophlan.tar.gz uniref.tar.gz utility.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")
check_file "humann: ChocoPhlAn pangenomes in chocophlan/" "humann_db/chocophlan/SGB1871_pangenome90.fna.gz"
check_file "humann: DIAMOND database in uniref/" "humann_db/uniref/demo_proteins_uniref90_uniref50_like_v2019_06.dmnd"
check_file "humann: KO map in utility_mapping/" "humann_db/utility_mapping/map_ko_uniref90.txt.gz"
check_file "humann: MetaPhlAn .pkl in metaphlan/" "humann_db/metaphlan/${IDX}.pkl"
check_file "humann: Bowtie2 index brought up from the tar's folder" "humann_db/metaphlan/${IDX}.rev.2.bt2l"
if [ "$(cat "${LAST_WORKDIR}/humann_db/DB_VERSION" 2>/dev/null)" = "v4_alpha-test: HUMAnN-test" ]; then
    echo "✓ humann: single-line DB_VERSION records the registry key (ChocoPhlAn flavour) and release"
    PASS=$((PASS + 1))
else
    echo "✗ humann: DB_VERSION missing or wrong"
    FAIL=$((FAIL + 1))
fi
if [ -z "$(ls "${LAST_WORKDIR}" | grep -v '^humann_db$')" ]; then
    echo "✓ humann: downloaded archives removed after unpacking"
    PASS=$((PASS + 1))
else
    echo "✗ humann: leftover files next to humann_db/"
    FAIL=$((FAIL + 1))
fi
expect_pass "humann: ChocoPhlAn archive wrapped in a folder is flattened" humann_db_reformat.sh \
    $(hu_args chocophlan_wrapped.tar.gz uniref.tar.gz utility.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")
check_file "humann: wrapped ChocoPhlAn flattened into chocophlan/" "humann_db/chocophlan/SGB2091_pangenome90.fna.gz"
expect_fail "humann: Bowtie2 index md5 mismatch fails" humann_db_reformat.sh \
    $(hu_args chocophlan.tar.gz uniref.tar.gz utility.tar.gz "${IDX}_bt2.tar" bad_bt2.md5)
expect_fail "humann: incomplete Bowtie2 index fails" humann_db_reformat.sh \
    $(hu_args chocophlan.tar.gz uniref.tar.gz utility.tar.gz bt2_incomplete.tar bt2_incomplete.md5)
expect_fail "humann: utility mapping without the KO map fails" humann_db_reformat.sh \
    $(hu_args chocophlan.tar.gz uniref.tar.gz utility_noko.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")
expect_fail "humann: ChocoPhlAn archive without pangenomes fails" humann_db_reformat.sh \
    $(hu_args empty.tar.gz uniref.tar.gz utility.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")
expect_fail "humann: UniRef archive without a .dmnd fails" humann_db_reformat.sh \
    $(hu_args chocophlan.tar.gz chocophlan.tar.gz utility.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")
expect_fail "humann: truncated archive fails" humann_db_reformat.sh \
    $(hu_args truncated.tar.gz uniref.tar.gz utility.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")
expect_fail "humann: 404 archive URL fails" humann_db_reformat.sh \
    $(hu_args no_such.tar.gz uniref.tar.gz utility.tar.gz "${IDX}_bt2.tar" "${IDX}_bt2.md5")

#
# Report image smoke test (audit #21): all libraries importable, no runtime pip
#
echo "--- report container image ---"
if docker run --rm "${REPORT_IMG}" python3 -c "import pandas, numpy, matplotlib, seaborn, h5py, biom" > /dev/null 2>&1; then
    echo "✓ report image: pandas/numpy/matplotlib/seaborn/h5py/biom importable"
    PASS=$((PASS + 1))
else
    echo "✗ report image: python imports failed"
    FAIL=$((FAIL + 1))
fi

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
