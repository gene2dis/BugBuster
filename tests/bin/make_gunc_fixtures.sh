#!/bin/bash
#
# Regenerates tests/data/gunc/gunc_db_progenomes2.1.dmnd, the stand-in GUNC
# reference of tests/modules/gunc_run.nf.test: a DIAMOND database of the first
# 20 Prodigal proteins of tests/data/bins/real_bins/bin.1.fa, built with the
# DIAMOND of the pinned GUNC image (2.1.24, which GUNC requires). Its sequence
# names are no reference genome GUNC knows, so every bin is mapped but left
# unscored (nan): enough to check the module's names and layout. The file name
# is GUNC's, from which GUNC reads the reference (proGenomes 2.1).
#
set -euo pipefail
REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd )"
IMG=$(grep -o "'quay.io/biocontainers/gunc:[^']*'" "${REPO_DIR}/modules/local/gunc_run/main.nf" | tr -d "'")
docker run --rm -u "$(id -u):$(id -g)" -v "${REPO_DIR}:${REPO_DIR}" -w "${REPO_DIR}/tests/data/gunc" "${IMG}" sh -c '
    prodigal -i ../bins/real_bins/bin.1.fa -a /tmp/p.faa -o /dev/null -p meta -q
    awk "/^>/ { n++ } n <= 20" /tmp/p.faa | sed "s/\*//g" > /tmp/small.faa
    diamond makedb --in /tmp/small.faa --db gunc_db_progenomes2.1 --quiet'
