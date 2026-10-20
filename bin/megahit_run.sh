#!/usr/bin/env bash
#
# megahit_run.sh <memory bytes|-> <megahit arguments...>
#
# Runs MEGAHIT for the MEGAHIT process with a memory budget and reports
# running out of memory as exit 137, the status the pipeline's retry strategy
# reruns with more memory.
#
# - Without -m, MEGAHIT sizes its graph build (SdBG) to 90% of the host's RAM,
#   not the task's limit: on a large host the container kills it at that
#   limit. The first argument is the task memory in bytes; MEGAHIT gets 90% of
#   it (the rest of MEGAHIT works outside that budget). '-' passes no -m. A -m
#   among the MEGAHIT arguments (ext.args) comes later and wins.
# - -m is a budget, not a limit: when a step needs more than -m, MEGAHIT
#   1.2.9 uses that step's minimum anyway (logged as "Minimum memory
#   required"), and the container kills megahit_core if it exceeds the task
#   memory. MEGAHIT relays that -9 as exit 247, which is not in the retry
#   list, so the task failed at once instead of getting more memory; 247 and
#   megahit_core's own fatal "No enough memory" (exit 1, in <out dir>/log)
#   now exit 137. Other failures keep their status.
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: megahit_run.sh <memory bytes|-> <megahit arguments...>" >&2
    exit 2
fi

mem_bytes="$1"
shift

mem_args=()
if [ "${mem_bytes}" != "-" ]; then
    if ! [[ "${mem_bytes}" =~ ^[0-9]+$ ]]; then
        echo "megahit_run.sh: memory must be a number of bytes or '-', got '${mem_bytes}'" >&2
        exit 2
    fi
    mem_args=(-m $(( mem_bytes * 9 / 10 )))
fi

# The output directory (-o/--out-dir) holds MEGAHIT's log
out_dir="megahit_out"
args=("$@")
for (( i = 0; i < ${#args[@]}; i++ )); do
    case "${args[i]}" in
        -o|--out-dir) out_dir="${args[i + 1]:-${out_dir}}" ;;
        --out-dir=*)  out_dir="${args[i]#--out-dir=}" ;;
    esac
done

rc=0
megahit "${mem_args[@]}" "$@" || rc=$?

if [ "${rc}" -ne 0 ]; then
    if [ "${rc}" -eq 247 ] || grep -qs "No enough memory" "${out_dir}/log"; then
        echo "megahit_run.sh: MEGAHIT ran out of memory (exit status ${rc}); exiting 137 so the task is retried with more memory" >&2
        exit 137
    fi
    exit "${rc}"
fi
