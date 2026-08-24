#!/bin/bash
#
# identify_bin_dir.sh <bin_dir_path> <metawrap_dir_suffix>
#
# Recognizes a binner output directory by name and prints "<sample_id>\t<binner>".
# The MetaWRAP directory name depends on params.metawrap_completeness/contamination
# (e.g. "metawrap_50_10_bins"), so the expected suffix is passed in by the caller
# instead of being hardcoded here (audit #6).
#
# Exits 1 on an unrecognized directory name: the batch QC/taxonomy processes only
# receive binner output dirs, so an unknown name means a naming drift bug and must
# fail loudly rather than silently yield empty reports.
#
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: identify_bin_dir.sh <bin_dir_path> <metawrap_dir_suffix>" >&2
    exit 2
fi

bin_dir="$1"
metawrap_suffix="$2"
name="$(basename "$bin_dir")"

if [ -n "$metawrap_suffix" ] && [[ "$name" == *_"$metawrap_suffix" ]]; then
    printf '%s\t%s\n' "${name%_"$metawrap_suffix"}" "metawrap"
elif [[ "$name" == *_metabat_bins ]]; then
    printf '%s\t%s\n' "${name%_metabat_bins}" "metabat"
elif [[ "$name" == *_semibin_output_bins ]]; then
    printf '%s\t%s\n' "${name%_semibin_output_bins}" "semibin"
elif [[ "$name" == *_comebin_bins ]]; then
    printf '%s\t%s\n' "${name%_comebin_bins}" "comebin"
elif [[ "$name" == comebin_res_bins ]]; then
    # COMEBin emits the nested <id>_comebin_bins/comebin_res/comebin_res_bins dir,
    # so the staged basename carries no sample id — recover it from the resolved
    # path, which still contains the <id>_comebin_bins component.
    target="$(readlink -f "$bin_dir")"
    if [[ "$target" == *_comebin* ]]; then
        prefix="${target%%_comebin*}"
        printf '%s\t%s\n' "${prefix##*/}" "comebin"
    else
        echo "ERROR: cannot derive sample id for COMEBin dir '$bin_dir' (resolved: '$target')" >&2
        exit 1
    fi
else
    echo "ERROR: unrecognized bin directory name '$name' (expected *_metabat_bins, *_semibin_output_bins, *_comebin_bins, comebin_res_bins, or *_${metawrap_suffix})" >&2
    exit 1
fi
