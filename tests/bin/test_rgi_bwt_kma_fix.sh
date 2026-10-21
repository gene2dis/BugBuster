#!/bin/bash
#
# Regression tests for bin/rgi_bwt_kma_fix.py (RGI_BWT). RGI's
# BWT.align_kma writes every read KMA reports, unplaced ones included, to the
# SAM that all later rgi bwt steps read (~390 bytes per clean read on disk,
# and 6.0.8 loads the read table into memory). Runs inside the pinned RGI
# image (the patch imports RGI's own module) with its real KMA and checks:
#   - the version guard accepts the installed 6.0.8 methods and refuses any
#     other one;
#   - on paired reads where one pair maps and one does not, the unpatched
#     method writes the unmapped reads to the SAM and the patched one writes
#     only the mapped reads, with the same alignments;
#   - the same for the interleaved (single-file) method;
#   - the filter drops records KMA writes at POS 0 (k-mer hits without an
#     alignment, unplaced once samtools reads them) and keeps placed mates;
#   - the kept records are grouped by read name (mates on adjacent lines, for
#     samtools fixmate), keeping their order within a name;
#   - a failing KMA stops the patched method; the unpatched one ignores it.
#
# Requirements: docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
RGI_IMG="quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

PASS=0
FAIL=0

check() {
    # check <desc> <command...>
    local desc="$1"
    shift
    if "$@"; then
        echo "✓ ${desc}"
        PASS=$((PASS + 1))
    else
        echo "✗ ${desc}"
        FAIL=$((FAIL + 1))
    fi
}

in_rgi() {
    docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -e MPLCONFIGDIR=/tmp \
        -v "${REPO_DIR}:${REPO_DIR}:ro" -v "${TMP_DIR}:${TMP_DIR}" -w "${TMP_DIR}" \
        "${RGI_IMG}" "$@"
}

echo "=== rgi_bwt_kma_fix.py test suite ==="

cat > "${TMP_DIR}/run_tests.py" <<'PY'
import importlib.util
import os
import random
import subprocess
import sys

shim_path, work = sys.argv[1], sys.argv[2]
os.chdir(work)
spec = importlib.util.spec_from_file_location('rgi_bwt_kma_fix', shim_path)
shim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shim)

import app.BWT
from app.BWT import BWT
shim.logger = app.BWT.logger

results = []
def check(desc, ok):
    results.append(ok)
    print(('PASS ' if ok else 'FAIL ') + desc)

# A 1 kb reference indexed with KMA; pair_hit is read from it, pair_miss is
# random sequence that maps nowhere
random.seed(1)
rand = lambda n: ''.join(random.choice('ACGT') for _ in range(n))
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
revcomp = lambda s: ''.join(comp[c] for c in reversed(s))
ref = rand(1000)
with open('ref.fasta', 'w') as fh:
    fh.write('>ARO:3000001|ID:1|Name:geneA\n' + ref + '\n')
subprocess.run(['kma', 'index', '-i', 'ref.fasta', '-o', 'idx'], check=True,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
pairs = {'pair_hit': (ref[100:250], revcomp(ref[400:550])), 'pair_miss': (rand(150), rand(150))}
for mate, path in ((0, 'r1.fastq'), (1, 'r2.fastq')):
    with open(path, 'w') as fh:
        for name, seqs in pairs.items():
            fh.write(f'@{name}/{mate + 1}\n{seqs[mate]}\n+\n{"I" * 150}\n')
with open('int.fastq', 'w') as fh:
    for name, seqs in pairs.items():
        for mate in (0, 1):
            fh.write(f'@{name}/{mate + 1}\n{seqs[mate]}\n+\n{"I" * 150}\n')

def bwt(read_one, read_two=None):
    obj = object.__new__(BWT)
    obj.read_one, obj.read_two, obj.threads = read_one, read_two, 1
    obj.check_index = lambda index_directory, reference_genome: None
    return obj

def records(sam):
    with open(sam) as fh:
        return [l.split('\t') for l in fh if not l.startswith('@')]

def names(recs):
    return sorted(r[0].split('/')[0] for r in recs)

# Guard: refuses changed methods, then patches the installed ones
originals = {n: getattr(BWT, n) for n in ('align_kma', 'align_kma_interleaved')}
class Changed(object):
    def align_kma(self):
        return None
    def align_kma_interleaved(self):
        return None
try:
    shim.patch(Changed)
    check('guard refuses a changed method', False)
except SystemExit as e:
    check('guard refuses a changed method', 'not the 6.0.8 method' in str(e))

# Unpatched 6.0.8 methods (the bug): unmapped reads go to the SAM
originals['align_kma'](bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'idx', 'orig.sam')
orig = records('orig.sam')
check('unpatched align_kma writes the unmapped pair to the SAM',
      'pair_miss' in names(orig) and any(int(r[1]) & 4 for r in orig))
originals['align_kma_interleaved'](bwt('int.fastq'), 'ref.fasta', 'idx', 'orig_int.sam')
check('unpatched align_kma_interleaved writes the unmapped pair to the SAM',
      'pair_miss' in names(records('orig_int.sam')))

shim.patch(BWT)
check('guard accepts the installed RGI 6.0.8 methods',
      BWT.align_kma is shim.align_kma and BWT.align_kma_interleaved is shim.align_kma_interleaved)

# Patched: mapped reads only, same alignments as the unpatched SAM
BWT.align_kma(bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'idx', 'fixed.sam')
fixed = records('fixed.sam')
check('patched align_kma writes only mapped reads',
      names(fixed) == ['pair_hit', 'pair_hit'] and not any(int(r[1]) & 4 for r in fixed))
check('patched align_kma keeps the unpatched alignments of the mapped reads',
      sorted(map(tuple, fixed)) == sorted(tuple(r) for r in orig if not int(r[1]) & 4))
BWT.align_kma_interleaved(bwt('int.fastq'), 'ref.fasta', 'idx', 'fixed_int.sam')
check('patched align_kma_interleaved writes only mapped reads',
      names(records('fixed_int.sam')) == ['pair_hit', 'pair_hit'])

# The filter: KMA writes reads with k-mer hits but no alignment at POS 0 with
# a reference name; samtools reads them as unplaced unmapped reads
ref_name = 'ARO:3000001|ID:1|Name:geneA'
sam = ['@HD\tVN:1.6', f'@SQ\tSN:{ref_name}\tLN:1000',
       f'placed/1\t99\t{ref_name}\t101\t60\t150M\t=\t401\t450\tA\tI',
       f'unmapped_mate/2\t165\t{ref_name}\t101\t0\t*\t=\t101\t0\tA\tI',
       f'pos0/1\t69\t{ref_name}\t0\t0\t*\t*\t0\t0\tA\tI',
       'nohit/1\t77\t*\t0\t0\t*\t*\t0\t0\tA\tI']
kept = subprocess.run(['awk', '-F', '\t', shim.PLACED_READS_AWK], input='\n'.join(sam) + '\n',
                      capture_output=True, text=True, check=True).stdout.splitlines()
check('filter keeps the header, placed reads and placed unmapped mates; drops POS 0 and unplaced reads',
      [l.split('\t')[0] for l in kept] == ['@HD', '@SQ', 'placed/1', 'unmapped_mate/2'])

# Grouping: KMA writes records in thread order, mates often apart; fixmate
# pairs adjacent records only
with open('group.sam', 'w') as fh:
    fh.write('\n'.join(['@HD\tVN:1.6', f'@SQ\tSN:{ref_name}\tLN:1000',
                        f'b\t99\t{ref_name}\t5\t60\t1M\t*\t0\t0\tA\tI',
                        f'a\t83\t{ref_name}\t9\t60\t1M\t*\t0\t0\tA\tI',
                        f'b\t147\t{ref_name}\t7\t60\t1M\t*\t0\t0\tA\tI',
                        f'a\t163\t{ref_name}\t3\t60\t1M\t*\t0\t0\tA\tI']) + '\n')
shim.group_by_name('group.sam')
with open('group.sam') as fh:
    grouped = [l.split('\t')[:2] for l in fh.read().splitlines()]
check('kept records grouped by read name, header first, order kept within a name',
      grouped == [['@HD', 'VN:1.6'], ['@SQ', f'SN:{ref_name}'],
                  ['a', '83'], ['a', '163'], ['b', '99'], ['b', '147']])

# A failing KMA (missing index): the unpatched method returns, the patched one stops
originals['align_kma'](bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'no_such_index', 'bad_orig.sam')
check('unpatched align_kma ignores a failing KMA', True)
try:
    BWT.align_kma(bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'no_such_index', 'bad.sam')
    check('patched align_kma stops when KMA fails', False)
except SystemExit as e:
    check('patched align_kma stops when KMA fails', 'KMA failed' in str(e))

sys.exit(0 if all(results) else 1)
PY

while IFS= read -r line; do
    case "${line}" in
        PASS\ *) echo "✓ ${line#PASS }"; PASS=$((PASS + 1)) ;;
        FAIL\ *) echo "✗ ${line#FAIL }"; FAIL=$((FAIL + 1)) ;;
    esac
done < <(in_rgi python3 -B run_tests.py "${REPO_DIR}/bin/rgi_bwt_kma_fix.py" "${TMP_DIR}" 2>/dev/null || true)

wrapper_help() {
    in_rgi python3 -B "${REPO_DIR}/bin/rgi_bwt_kma_fix.py" bwt -h | grep -q 'usage: rgi bwt'
}
check "wrapper patches and hands over to rgi bwt (--help)" wrapper_help

echo ""
echo "=== Results: ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ] && [ "${PASS}" -ge 12 ]
