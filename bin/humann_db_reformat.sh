#!/bin/bash
#
# Provision the HUMAnN 4.0.0a2 database set for the read-level functional
# branch (design doc Section 4.6.1, task T8c, Q2).
# Usage: humann_db_reformat.sh <key> <chocophlan-url> <uniref-url> <utility-url>
#            <metaphlan-index> <metaphlan-url> <metaphlan-md5-url>
#            <metaphlan-bt2-url> <metaphlan-bt2-md5-url> <provenance>
#   key                 registry key (v4_alpha-full | v4_alpha-ec_filtered),
#                       recorded in DB_VERSION (the ChocoPhlAn flavour)
#   chocophlan-url      ChocoPhlAn v4_alpha tar.gz (full or EC-filtered)
#   uniref-url          UniRef90 EC-filtered DIAMOND tar.gz (only v4 protein DB)
#   utility-url         full_mapping_v4_alpha tar.gz (KO/EC maps, MetaCyc files)
#   metaphlan-index     the MetaPhlAn database tag HUMAnN 4.0.0a2 accepts
#                       (mpa_vOct22_CHOCOPhlAnSGB_202403)
#   metaphlan-url       its base tar (.pkl + marker files)
#   metaphlan-md5-url   the published md5 of that tar
#   metaphlan-bt2-url   its prebuilt Bowtie2 (.bt2l) index tar
#   metaphlan-bt2-md5-url  the published md5 of that tar
#   provenance          human-readable release string for DB_VERSION
#
# Produces humann_db/, the database ROOT the HUMANN module expects:
#   humann_db/chocophlan/       nucleotide pangenome database (--nucleotide-database)
#   humann_db/uniref/           *.dmnd translated-search database (--protein-database)
#   humann_db/utility_mapping/  map_*.txt.gz + the MetaCyc pathway files (--utility-database,
#                               --pathways-database)
#   humann_db/metaphlan/        <index>.pkl + <index>.*.bt2l (METAPHLAN_DB_DIR)
#   humann_db/DB_VERSION        single line (read by HUMANN into versions.yml)
# HUMAnN publishes no checksums for its archives, so those get non-empty and
# layout checks; the MetaPhlAn tars are md5-verified. Archives are deleted
# after extraction. Each archive is extracted into its own folder (as
# humann_databases does); a single wrapping subdirectory is flattened.
#
set -euo pipefail

usage="usage: humann_db_reformat.sh <key> <chocophlan-url> <uniref-url> <utility-url> <metaphlan-index> <metaphlan-url> <metaphlan-md5-url> <metaphlan-bt2-url> <metaphlan-bt2-md5-url> <provenance>"
key="${1:?${usage}}"
chocophlan_url="${2:?${usage}}"
uniref_url="${3:?${usage}}"
utility_url="${4:?${usage}}"
mpa_index="${5:?${usage}}"
mpa_url="${6:?${usage}}"
mpa_md5_url="${7:?${usage}}"
bt2_url="${8:?${usage}}"
bt2_md5_url="${9:?${usage}}"
provenance="${10:?${usage}}"

# Overridable for the test harness (local HTTP server, no waits)
retries="${HUMANN_DB_RETRIES:-5}"
retry_wait="${HUMANN_DB_RETRY_WAIT:-10}"

fetch() {  # fetch <url> <output>
    wget --tries="${retries}" --continue --waitretry="${retry_wait}" --no-verbose -O "$2" "$1"
}

expected_md5() {  # expected_md5 <md5-url> <tar-name> -> prints the md5 for that file
    local listing
    listing="$(wget --tries="${retries}" --waitretry="${retry_wait}" -q -O - "$1")"
    awk -v name="$2" '$2 == name || $2 == "*" name { print $1; exit }' <<< "${listing}"
}

verify_md5() {  # verify_md5 <file> <expected> <label>
    local actual
    actual="$(md5sum "$1" | cut -d' ' -f1)"
    if [ -z "$2" ]; then
        echo "ERROR: no published md5 found for $3" >&2
        exit 1
    fi
    if [ "${actual}" != "$2" ]; then
        echo "ERROR: md5 mismatch for $3 (got ${actual}, expected $2) - corrupt or truncated download" >&2
        exit 1
    fi
}

flatten() {  # flatten <dir>: if <dir> holds only one subdirectory, move its contents up
    local dir="$1" entries
    entries=("${dir}"/*)
    if [ "${#entries[@]}" -eq 1 ] && [ -d "${entries[0]}" ]; then
        local inner="${entries[0]}"
        find "${inner}" -mindepth 1 -maxdepth 1 -exec mv -t "${dir}" {} +
        rmdir "${inner}"
    fi
}

extract_into() {  # extract_into <url> <dir> <label>
    local archive="${2}.download"
    mkdir -p "$2"
    fetch "$1" "${archive}"
    if ! tar -xf "${archive}" -C "$2"; then
        echo "ERROR: could not extract the $3 archive ($1) - corrupt or truncated download" >&2
        exit 1
    fi
    rm -f "${archive}"
    flatten "$2"
}

root=humann_db
mkdir -p "${root}"

# --- HUMAnN databases (no published checksums) ---
extract_into "${chocophlan_url}" "${root}/chocophlan" ChocoPhlAn
if [ -z "$(find "${root}/chocophlan" -type f -name '*.gz' -size +0 -print -quit)" ]; then
    echo "ERROR: ${root}/chocophlan holds no non-empty pangenome (*.gz) files - not a ChocoPhlAn archive" >&2
    exit 1
fi

extract_into "${uniref_url}" "${root}/uniref" UniRef90
if [ -z "$(find "${root}/uniref" -type f -name '*.dmnd' -size +0 -print -quit)" ]; then
    echo "ERROR: ${root}/uniref holds no non-empty DIAMOND (*.dmnd) database - not a UniRef archive" >&2
    exit 1
fi

extract_into "${utility_url}" "${root}/utility_mapping" "utility mapping"
# the files the pipeline uses (KO/EC regrouping + names, MetaCyc pathways)
for required in map_ko_uniref90.txt.gz map_level4ec_uniclust90.txt.gz map_ko_name.txt.gz \
                map_level4ec_name.txt.gz metacyc_reactions_level4ec_only.uniref.bz2 \
                metacyc_pathways_structured_filtered_v24_subreactions; do
    if [ ! -s "${root}/utility_mapping/${required}" ]; then
        echo "ERROR: ${root}/utility_mapping/${required} missing or empty - not the full_mapping_v4_alpha archive" >&2
        exit 1
    fi
done

# --- MetaPhlAn database (md5-verified) ---
mkdir -p "${root}/metaphlan"
mpa_tar="${mpa_index}.tar"
bt2_tar="${mpa_index}_bt2.tar"
fetch "${mpa_url}" "${mpa_tar}"
verify_md5 "${mpa_tar}" "$(expected_md5 "${mpa_md5_url}" "${mpa_tar}")" "${mpa_tar}"
fetch "${bt2_url}" "${bt2_tar}"
verify_md5 "${bt2_tar}" "$(expected_md5 "${bt2_md5_url}" "${bt2_tar}")" "${bt2_tar}"
for archive in "${mpa_tar}" "${bt2_tar}"; do
    tar -xf "${archive}" -C "${root}/metaphlan"
    rm -f "${archive}"
done
flatten "${root}/metaphlan"
# the two tars may unpack into different subfolders: bring every file of the
# index to the top level, where MetaPhlAn (METAPHLAN_DB_DIR) looks for it
find "${root}/metaphlan" -mindepth 2 -type f -name "${mpa_index}*" -exec mv -t "${root}/metaphlan" {} +
if [ ! -s "${root}/metaphlan/${mpa_index}.pkl" ]; then
    echo "ERROR: ${root}/metaphlan/${mpa_index}.pkl missing or empty after unpacking" >&2
    exit 1
fi
for ext in 1.bt2l 2.bt2l 3.bt2l 4.bt2l rev.1.bt2l rev.2.bt2l; do
    if [ ! -s "${root}/metaphlan/${mpa_index}.${ext}" ]; then
        echo "ERROR: ${root}/metaphlan/${mpa_index}.${ext} missing or empty - incomplete MetaPhlAn Bowtie2 index" >&2
        exit 1
    fi
done

printf '%s\n' "${key}: ${provenance}" > "${root}/DB_VERSION"
