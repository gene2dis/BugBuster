#!/bin/bash
#
# Download the Web of Life release 2 (WoLr2) files needed by the Woltka
# read-level functional backend (design doc Section 4.6.2, task T8a).
# Usage: woltka_db_reformat.sh <wol2-base-url>   (e.g. https://ftp.microbio.me/pub/wol2)
#
# Produces woltka_db/ mirroring the FTP layout verbatim (only the files the
# backend reads), so a hand-made mirror of the same files is a valid
# --custom_woltka_db:
#   databases/bowtie2/WoLr2.{1,2,3,4,rev.1,rev.2}.bt2l   Bowtie2 index, 93.6 GB
#   proteins/{coords.txt.xz,length.map.xz}               ORF coordinates/lengths
#   function/kegg/{orf-to-ko.map.xz,ko-to-ec.map,ko-to-cog.map,ko_name.txt}
#   function/metacyc/{orf-to-protein.map.xz,protein-to-enzrxn.map,
#       enzrxn-to-reaction.map,reaction-to-pathway.map,pathway_name.txt}
#   function/pfam/{orf-to-pfam.map.xz,pfam_name.txt}
# (~94 GB total). The pipeline bundles none of this data: KEGG- and
# MetaCyc-derived maps are fetched from the official public host on the
# user's behalf (owner decision, design doc Q5).
#
# Integrity: every .xz file that WoLr2 publishes an .md5 for is verified
# against it - those md5s are of the UNCOMPRESSED content (verified
# 2026-10-07), hence `xz -dc | md5sum`. The Bowtie2 index ships no md5; its
# files are checked non-empty (a corrupt index fails loudly at alignment).
# Runs in the pinned woltka container (FORMAT_WOLTKA_DB), whose wget is
# busybox (no --tries/--waitretry), so retries are looped here; -c resumes a
# partial file. A DB_VERSION file records the release for provenance
# (WOLTKA_CLASSIFY reads it into versions.yml).
#
set -euo pipefail

base_url="${1:?usage: woltka_db_reformat.sh <wol2-base-url>}"
base_url="${base_url%/}"

index_files="WoLr2.1.bt2l WoLr2.2.bt2l WoLr2.3.bt2l WoLr2.4.bt2l WoLr2.rev.1.bt2l WoLr2.rev.2.bt2l"
# files with a published md5 of their uncompressed content
md5_files="proteins/coords.txt proteins/length.map function/kegg/orf-to-ko.map function/metacyc/orf-to-protein.map function/pfam/orf-to-pfam.map"
plain_files="function/kegg/ko-to-ec.map function/kegg/ko-to-cog.map function/kegg/ko_name.txt function/metacyc/protein-to-enzrxn.map function/metacyc/enzrxn-to-reaction.map function/metacyc/reaction-to-pathway.map function/metacyc/pathway_name.txt function/pfam/pfam_name.txt"

# Retry policy (overridable for the test harness, which serves fixtures
# locally and must not wait on deliberate 404s)
retries="${WOLTKA_DB_RETRIES:-5}"
retry_wait="${WOLTKA_DB_RETRY_WAIT:-10}"

fetch() {
    local rel="$1" attempt
    mkdir -p "woltka_db/$(dirname "${rel}")"
    for attempt in $(seq 1 "${retries}"); do
        if wget -c -q -O "woltka_db/${rel}" "${base_url}/${rel}"; then
            return 0
        fi
        echo "WARN: download of ${rel} failed (attempt ${attempt}/${retries}), retrying" >&2
        sleep "${retry_wait}"
    done
    echo "ERROR: could not download ${base_url}/${rel}" >&2
    exit 1
}

for f in ${index_files}; do
    fetch "databases/bowtie2/${f}"
done
for f in ${md5_files}; do
    fetch "${f}.xz"
    fetch "${f}.md5"
done
for f in ${plain_files}; do
    fetch "${f}"
done

# Verify before finishing, so the exit status reflects the database
for f in ${index_files}; do
    if [ ! -s "woltka_db/databases/bowtie2/${f}" ]; then
        echo "ERROR: expected index file 'woltka_db/databases/bowtie2/${f}' missing or empty after download" >&2
        exit 1
    fi
done
for f in ${md5_files}; do
    expected="$(cut -d' ' -f1 "woltka_db/${f}.md5")"
    actual="$(xz -dc "woltka_db/${f}.xz" | md5sum | cut -d' ' -f1)"
    if [ -z "${expected}" ] || [ "${expected}" != "${actual}" ]; then
        echo "ERROR: md5 mismatch for woltka_db/${f}.xz (uncompressed md5 ${actual}, published ${expected:-missing}) - corrupt or truncated download" >&2
        exit 1
    fi
done
for f in ${plain_files}; do
    if [ ! -s "woltka_db/${f}" ]; then
        echo "ERROR: expected file 'woltka_db/${f}' missing or empty after download" >&2
        exit 1
    fi
done

printf '%s\n' "WoLr2 (${base_url})" > woltka_db/DB_VERSION
