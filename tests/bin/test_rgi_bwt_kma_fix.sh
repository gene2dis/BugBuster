#!/bin/bash
#
# Regression tests for bin/rgi_bwt_kma_fix.py (RGI_BWT). RGI's
# BWT.align_* methods write every read the aligner reports, unplaced ones
# included, to the SAM that all later rgi bwt steps read (~390 bytes per clean
# read on disk, and 6.0.8 loads the read table into memory). Runs inside the
# pinned RGI image (the patch imports RGI's own module) with its real KMA,
# bowtie2 and bwa and checks:
#   - the version guard accepts the installed 6.0.8 methods and refuses any
#     other one;
#   - on paired reads where one pair maps and one does not, the unpatched
#     method writes the unmapped reads to the SAM and the patched one writes
#     only the mapped reads, with the same alignments;
#   - the same for the interleaved (single-file) method;
#   - the filter drops records KMA writes at POS 0 (k-mer hits without an
#     alignment, unplaced once samtools reads them) and keeps placed mates;
#   - the kept records are grouped by read name (mates on adjacent lines, for
#     samtools fixmate), keeping their order within a name, and the @HD line
#     says so; a SAM without placed records keeps its header;
#   - a failing KMA stops the patched method; the unpatched one ignores it;
#   - the bowtie2 and bwa copies drop the pair that maps nowhere, keep every
#     other record of the unpatched SAM, and stop when the aligner fails;
#   - a killed aligner or grouping sort stops rgi with 128 + signal (137, so
#     RGI_BWT retries with more memory), and a failed filter is reported as
#     the filter, not as the aligner it took down with SIGPIPE;
#   - rebuilding an incomplete index replaces the task's symlinks to it and
#     leaves the linked (shared database) files untouched.
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
import contextlib
import importlib.util
import io
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
for index_cmd in (['kma', 'index', '-i', 'ref.fasta', '-o', 'idx'],
                  ['bowtie2-build', '--quiet', 'ref.fasta', 'bt2idx'],
                  ['bwa', 'index', '-p', 'bwaidx', 'ref.fasta']):
    subprocess.run(index_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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

def stops(call):
    """Run call; return (exit code, stderr) if it stops rgi, else None."""
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err):
            call()
    except SystemExit as e:
        return e.code, err.getvalue()
    return None

# Guard: refuses changed methods, then patches the installed ones
originals = {n: getattr(BWT, n) for n, _ in shim.REPLACEMENTS}
class Changed(object):
    pass
for n in originals:
    setattr(Changed, n, lambda self: None)
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
      all(getattr(BWT, n) is f for n, f in shim.REPLACEMENTS)
      and BWT.create_index.__qualname__ == 'patch.<locals>.create_index')

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
def kma_filter(lines):
    """Run a SAM through run_kma (a fake KMA that prints it) and return the kept lines."""
    with open('in.sam', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    open('out.sam.temp.res', 'w').close()
    shim.run_kma('cat in.sam', 'out.sam')
    with open('out.sam') as fh:
        return fh.read().splitlines()
kept = kma_filter(sam)
check('filter keeps the header, placed reads and placed unmapped mates; drops POS 0 and unplaced reads',
      [l.split('\t')[0] for l in kept] == ['@HD', '@SQ', 'placed/1', 'unmapped_mate/2'])

# Grouping: KMA writes records in thread order, mates often apart; fixmate
# pairs adjacent records only
grouped = [l.split('\t')[:2] for l in kma_filter([
    '@HD\tVN:1.6\tGO:reference', f'@SQ\tSN:{ref_name}\tLN:1000',
    f'b\t99\t{ref_name}\t5\t60\t1M\t*\t0\t0\tA\tI',
    f'a\t83\t{ref_name}\t9\t60\t1M\t*\t0\t0\tA\tI',
    f'b\t147\t{ref_name}\t7\t60\t1M\t*\t0\t0\tA\tI',
    f'a\t163\t{ref_name}\t3\t60\t1M\t*\t0\t0\tA\tI'])]
check('kept records grouped by read name, header first, order kept within a name',
      grouped == [['@HD', 'VN:1.6'], ['@SQ', f'SN:{ref_name}'],
                  ['a', '83'], ['a', '163'], ['b', '99'], ['b', '147']])
check('the @HD line says the records are grouped by query name',
      kma_filter(['@HD\tVN:1.6\tGO:reference', 'a\t0\tr\t1\t60\t1M\t*\t0\t0\tA\tI'])[0] == '@HD\tVN:1.6\tGO:query')
check('a SAM without placed records keeps its header',
      kma_filter(['@HD\tVN:1.6', f'@SQ\tSN:{ref_name}\tLN:1000', 'nohit/1\t77\t*\t0\t0\t*\t*\t0\t0\tA\tI'])
      == ['@HD\tVN:1.6', f'@SQ\tSN:{ref_name}\tLN:1000'])

# bowtie2 and bwa: every record with a reference name, the rest unchanged
for method, index in (('align_bowtie2', 'bt2idx'), ('align_bwa_paired_end_mapping', 'bwaidx')):
    originals[method](bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', index, f'orig_{method}.sam')
    getattr(BWT, method)(bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', index, f'fixed_{method}.sam')
    orig, fixed = records(f'orig_{method}.sam'), records(f'fixed_{method}.sam')
    check(f'patched {method} drops the pair that maps nowhere and keeps the other records',
          'pair_miss' in names(orig) and names(fixed) == ['pair_hit', 'pair_hit']
          and fixed == [r for r in orig if r[2] != '*'])
for method in ('align_bowtie2', 'align_bwa_paired_end_mapping'):
    stop = stops(lambda: getattr(BWT, method)(
        bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'no_such_index', f'bad_{method}.sam'))
    check(f'patched {method} stops when the aligner fails',
          stop is not None and stop[0] != 0 and 'failed (exit status' in stop[1])

# A failing KMA (missing index): the unpatched method returns, the patched one stops
originals['align_kma'](bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'no_such_index', 'bad_orig.sam')
check('unpatched align_kma ignores a failing KMA', True)
stop = stops(lambda: BWT.align_kma(bwt('r1.fastq', 'r2.fastq'), 'ref.fasta', 'no_such_index', 'bad.sam'))
check('patched align_kma stops when KMA fails',
      stop is not None and stop[0] != 0 and 'KMA failed' in stop[1])

# Killed children: 128 + signal, the status RGI_BWT's errorStrategy retries
open('killed.sam.temp.res', 'w').close()
stop = stops(lambda: shim.run_kma('sh -c "kill -9 \\$\\$"', 'killed.sam'))
check('a killed KMA stops rgi with exit status 137',
      stop is not None and stop[0] == 137 and 'KMA failed' in stop[1])
os.makedirs('fake_sort', exist_ok=True)
for name, body in (('killed', 'cat >/dev/null; kill -9 $$'), ('early', 'exit 3')):
    os.makedirs(f'fake_sort/{name}', exist_ok=True)
    with open(f'fake_sort/{name}/sort', 'w') as fh:
        fh.write(f'#!/bin/sh\n{body}\n')
    os.chmod(f'fake_sort/{name}/sort', 0o755)
path = os.environ['PATH']
os.environ['PATH'] = f'{work}/fake_sort/killed:{path}'
stop = stops(lambda: kma_filter(sam))
check('a killed grouping sort stops rgi with exit status 137',
      stop is not None and stop[0] == 137 and 'filtering the KMA SAM failed' in stop[1])
# A sort that exits without reading: awk dies of SIGPIPE, then the aligner
with open('many.sam', 'w') as fh:
    fh.write('\n'.join(sam[:2] + [sam[2]] * 200000) + '\n')
os.environ['PATH'] = f'{work}/fake_sort/early:{path}'
open('early.sam.temp.res', 'w').close()
stop = stops(lambda: shim.run_kma('cat many.sam', 'early.sam'))
os.environ['PATH'] = path
check('a failed filter is reported as the filter, not as the aligner',
      stop is not None and stop[0] != 0 and 'filtering the KMA SAM failed' in stop[1])

# Rebuilding an incomplete staged index: the task's symlinks are replaced,
# the shared files they point to are not written
os.makedirs('shared/bwt', exist_ok=True)
os.makedirs('task/bwt', exist_ok=True)
with open('shared/bwt/kma.name', 'w') as fh:
    fh.write('shared\n')
os.symlink(os.path.abspath('shared/bwt/kma.name'), 'task/bwt/kma.name')
obj = bwt('r1.fastq', 'r2.fastq')
del obj.check_index
obj.aligner = 'kma'
BWT.create_index(obj, index_directory='task/bwt/kma', reference_genome='ref.fasta')
with open('shared/bwt/kma.name') as fh:
    shared_kept = fh.read() == 'shared\n'
check('an index rebuild writes new files in the task, not through the staged symlinks',
      shared_kept and not os.path.islink('task/bwt/kma.name')
      and os.path.isfile('task/bwt/kma.name') and os.path.isfile('task/bwt/kma.comp.b'))

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
[ "${FAIL}" -eq 0 ] && [ "${PASS}" -ge 22 ]
