#!/bin/bash
#
# Regression tests for audit items #3 and #10:
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

if [ -z "${DOWNLOAD_IMG}" ] || [ -z "${BAKTA_IMG}" ] || [ -z "${WOLTKA_IMG}" ] || [ -z "${REPORT_IMG}" ]; then
    echo "ERROR: could not extract pinned container images from the modules"
    exit 1
fi

echo "=== DB reformat script test suite ==="
echo "Download image: ${DOWNLOAD_IMG}"
echo "Bakta image:    ${BAKTA_IMG}"
echo "Woltka image:   ${WOLTKA_IMG}"
echo "Report image:   ${REPORT_IMG}"
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
        && tar -czf "${f}/gtdbtk_r232_data.tar.gz" release232 )
    mkdir -p "${staging}/gtdbtk_empty/release232"
    ( cd "${staging}/gtdbtk_empty" && tar -czf "${f}/gtdbtk_empty.tar.gz" release232 )

    # Mock NCBI blast db repository: metadata JSON + 2 volumes + md5 files
    mkdir -p "${f}/blastdb" "${staging}/blast"
    ( cd "${staging}/blast" \
        && echo t > taxdb.btd && echo t > taxdb.bti && echo t > taxonomy4blast.sqlite3 \
        && echo v > nt.000.nin && echo v > nt.000.nhr \
        && tar -czf "${f}/blastdb/nt.000.tar.gz" taxdb.btd taxdb.bti taxonomy4blast.sqlite3 nt.000.nin nt.000.nhr \
        && rm ./* \
        && echo v > nt.001.nin && echo v > nt.001.nhr \
        && tar -czf "${f}/blastdb/nt.001.tar.gz" nt.001.nin nt.001.nhr )
    ( cd "${f}/blastdb" \
        && md5sum nt.000.tar.gz > nt.000.tar.gz.md5 \
        && md5sum nt.001.tar.gz > nt.001.tar.gz.md5 \
        && printf '{"dbname":"nt","files":["%s","%s"]}\n' \
            "https://ftp.ncbi.nlm.nih.gov/blast/db/nt.000.tar.gz" \
            "https://ftp.ncbi.nlm.nih.gov/blast/db/nt.001.tar.gz" > nt-nucl-metadata.json )

    # Same blast repo but with a corrupted checksum for volume 001 ('x' can
    # never be a hex digit; the former '0' left the md5 unchanged whenever the
    # timestamp-dependent tarball md5 already started with 0 - a 1-in-16 flake)
    cp -r "${f}/blastdb" "${f}/blastdb_badmd5"
    sed 's/^./x/' "${f}/blastdb_badmd5/nt.001.tar.gz.md5" > "${f}/blastdb_badmd5/tmp.md5" \
        && mv "${f}/blastdb_badmd5/tmp.md5" "${f}/blastdb_badmd5/nt.001.tar.gz.md5"

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
# tax_files_reformat.sh
#
echo "--- tax_files_reformat.sh ---"
expect_pass "taxdump: valid tarball" tax_files_reformat.sh "${BASE_URL}/taxdump.tar.gz"
check_file "taxdump: nodes.dmp in tax_files" "tax_files/nodes.dmp"
check_file "taxdump: names.dmp in tax_files" "tax_files/names.dmp"
expect_fail "taxdump: corrupt tarball fails" tax_files_reformat.sh "${BASE_URL}/corrupt.tar.gz"
expect_fail "taxdump: tarball without dmp files fails" tax_files_reformat.sh "${BASE_URL}/k2_mini.tar.gz"

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
check_file "gtdbtk: release files extracted" "release232/markers/marker.txt"
expect_fail "gtdbtk: corrupt tarball fails" gtdb-tk_db_reformat.sh "${BASE_URL}/corrupt.tar.gz"
expect_fail "gtdbtk: empty package fails" gtdb-tk_db_reformat.sh "${BASE_URL}/gtdbtk_empty.tar.gz"

#
# blast_nt_reformat.sh
#
echo "--- blast_nt_reformat.sh ---"
expect_pass "blast: valid volume set" blast_nt_reformat.sh "${BASE_URL}/blastdb"
check_file "blast: taxdb.btd in blast_nt_db" "blast_nt_db/taxdb.btd"
check_file "blast: volume 000 index in blast_nt_db" "blast_nt_db/nt.000.nin"
check_file "blast: volume 001 index in blast_nt_db" "blast_nt_db/nt.001.nin"
expect_fail "blast: wrong md5 fails" blast_nt_reformat.sh "${BASE_URL}/blastdb_badmd5"
expect_fail "blast: missing metadata fails" blast_nt_reformat.sh "${BASE_URL}/no_such_dir"

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
