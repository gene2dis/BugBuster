#!/bin/bash
#
# bin/monitor_disk_usage.sh keeps running when du reports an error (roadmap
# Phase 21, RV12; review R6-06): a folder it cannot read (or a file removed
# while it scans, the normal state of a live work directory) made du exit
# non-zero while still printing its total, and the monitor stopped at its
# first iteration. Runs it for about three intervals on a work directory
# holding an unreadable folder and checks it is still running and logging.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
WORK=$(mktemp -d)
cleanup() { chmod -R u+rwx "${WORK}" 2>/dev/null || true; rm -rf "${WORK}"; }
trap cleanup EXIT
FAIL=0

mkdir -p "${WORK}/work/ok" "${WORK}/work/locked" "${WORK}/results"
head -c 4096 /dev/zero > "${WORK}/work/ok/file"
chmod 000 "${WORK}/work/locked"
if du -sb "${WORK}/work" > /dev/null 2>&1; then
    echo "SKIP: du can read an unreadable folder here (running as root?)"
    exit 0
fi

cd "${WORK}"
"${REPO_DIR}/bin/monitor_disk_usage.sh" "${WORK}/work" "${WORK}/results" 1 > monitor.out 2>&1 &
pid=$!
sleep 3.5
if kill -0 "${pid}" 2>/dev/null; then
    echo "✓ the monitor is still running after a du error"
    kill -TERM "${pid}" 2>/dev/null || true
    wait "${pid}" 2>/dev/null || true
else
    echo "✗ the monitor stopped: $(tail -n 3 monitor.out)"
    FAIL=1
fi
rows=$(($(wc -l < disk_usage_monitor.csv) - 1))
if [ "${rows}" -ge 2 ] && ! awk -F, 'NR > 1 && $2 !~ /^[0-9]+$/ { bad = 1 } END { exit !bad }' disk_usage_monitor.csv; then
    echo "✓ ${rows} CSV rows, each with one work-directory size"
else
    echo "✗ CSV rows: ${rows}"; cat disk_usage_monitor.csv
    FAIL=1
fi
[ "${FAIL}" -eq 0 ] && echo "PASS" || { echo "FAIL"; exit 1; }
