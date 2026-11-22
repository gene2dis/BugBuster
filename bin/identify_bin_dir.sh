#!/bin/bash
#
# identify_bin_dir.sh <bin_dir_path> [sample_id]
#
# Recognizes a binner output directory by name and prints "<sample_id>\t<binner>";
# Binette's refined bins (<sample>_binette_bins) report the binner "binette".
#
# The sample id is the second argument when the process passes it (meta.id,
# roadmap Phase 20, RV5): a named directory (<id>_metabat_bins, ...) must then
# carry that id, and COMEBin's nested comebin_res_bins directory, whose name
# holds no id, takes it from the argument alone. It used to be recovered from
# the resolved symlink target, which failed where inputs are copied (cloud
# executors) and misread a work path containing "_comebin".
#
# Exits 1 on an unrecognized directory name, on a name whose id differs from
# the given one, and on comebin_res_bins without an id: the batch QC/taxonomy
# processes only receive binner output dirs, so any of these is a naming or
# wiring bug and must fail loudly rather than silently yield empty reports.
#
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: identify_bin_dir.sh <bin_dir_path> [sample_id]" >&2
    exit 2
fi

bin_dir="$1"
given_id="${2:-}"
name="$(basename "$bin_dir")"

if [[ "$name" == *_binette_bins ]]; then
    named_id="${name%_binette_bins}"; binner="binette"
elif [[ "$name" == *_metabat_bins ]]; then
    named_id="${name%_metabat_bins}"; binner="metabat"
elif [[ "$name" == *_semibin_output_bins ]]; then
    named_id="${name%_semibin_output_bins}"; binner="semibin"
elif [[ "$name" == *_comebin_bins ]]; then
    named_id="${name%_comebin_bins}"; binner="comebin"
elif [[ "$name" == comebin_res_bins ]]; then
    # COMEBin emits <id>_comebin_bins/comebin_res/comebin_res_bins: the staged
    # name holds no sample id
    named_id=""; binner="comebin"
else
    echo "ERROR: unrecognized bin directory name '$name' (expected *_metabat_bins, *_semibin_output_bins, *_comebin_bins, comebin_res_bins or *_binette_bins)" >&2
    exit 1
fi

if [ -z "$named_id" ] && [ -z "$given_id" ]; then
    echo "ERROR: COMEBin directory '$bin_dir' carries no sample id: pass it as the second argument (meta.id)" >&2
    exit 1
fi
if [ -n "$named_id" ] && [ -n "$given_id" ] && [ "$named_id" != "$given_id" ]; then
    echo "ERROR: bin directory '$bin_dir' is named for sample '$named_id' but was staged for sample '$given_id'" >&2
    exit 1
fi

printf '%s\t%s\n' "${given_id:-$named_id}" "$binner"
