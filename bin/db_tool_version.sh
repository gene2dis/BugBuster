#!/bin/sh
#
# db_tool_version.sh <tool>
#
# The version of a tool a database preparation process (modules/local/format_db)
# runs, as one word for its versions.yml (roadmap Phase 29, CL1; issue #24).
# The images differ: GNU tools print "--version", BusyBox applets do not (the
# woltka image's wget, the busybox tar / gzip / xz / unzip of the biocontainers),
# so a BusyBox applet is reported as "busybox <version>". Never fails: the
# versions file is written before the download, and a version that cannot be
# read is "unknown", so it can never stop a run (or lose a download).
#
tool=${1:?usage: db_tool_version.sh <tool>}

command -v "$tool" > /dev/null 2>&1 || { echo "not found"; exit 0; }
case "$tool" in
    diamond)     out=$(diamond version 2>&1 | tr -d '\000') ;;
    makeblastdb) out=$(makeblastdb -version 2>&1 | tr -d '\000') ;;
    unzip)       out=$(unzip -v 2>&1 | tr -d '\000') ;;
    bzip2)       out=$(bzip2 --version 2>&1 < /dev/null | tr -d '\000') ;;
    *)           out=$($tool --version 2>&1 < /dev/null | tr -d '\000') ;;
esac
# a BusyBox applet (a symlink or a hard link to busybox) names BusyBox in its
# version or usage text ("tar (busybox) 1.36.1", "BusyBox v1.36.1 ...")
if printf '%s\n' "$out" | grep -qi busybox; then
    v=$( { busybox 2>&1 | head -n 1; printf '%s
' "$out"; } | grep -oE '[Vv]?[0-9]+(\.[0-9]+)+' | head -n 1 | tr -d vV)
    echo "busybox ${v:-unknown}"
    exit 0
fi
v=$(printf '%s\n' "$out" | head -n 1 | grep -oE '[0-9]+(\.[0-9]+)+[A-Za-z0-9+]*' | head -n 1)
[ -n "$v" ] || v=$(printf '%s\n' "$out" | grep -oE '[0-9]+(\.[0-9]+)+[A-Za-z0-9+]*' | head -n 1)
echo "${v:-unknown}"
