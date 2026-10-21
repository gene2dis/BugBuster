#!/bin/bash
#
# Regression tests for the count_bases patch in bin/run_microbe_census_py3fix.py
# (performance review M6). Runs inside the pinned MicrobeCensus image and checks:
#   - the version guard accepts the installed v1.1.1 code;
#   - the patched count_bases gives the stock total for plain 4-line FASTQ
#     (gzipped or not, last line without a newline, empty sequences, a
#     quality line starting with '@', an empty file) through its fast path;
#   - files the fast path must not read (CR line ends, FASTA, wrapped or
#     short-quality FASTQ, leading junk) fall back to the stock count;
#   - several files at once (the forked workers) sum to the stock total;
#   - the guard refuses a changed function.
#
# Requirements: docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
MC_IMG="quay.io/biocontainers/microbecensus:1.1.1--pyhca03a8a_2"

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

cd "${TMP_DIR}"
printf '@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n' > plain.fastq
gzip -c plain.fastq > plain.fastq.gz
printf '@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII' > no_last_newline.fastq
printf '@r1\n\n+\n\n@r2\nACG\n+\nIII\n' > empty_sequence.fastq
printf '@r1\nACGT\n+\n@III\n@r2\nAC\n+\nII\n' > quality_at.fastq
: > empty.fastq
printf '@r1\r\nACGT\r\n+\r\nIIII\r\n' > crlf.fastq
printf '>c1\nACGTACGT\nACG\n>c2\nAAA\n' > seqs.fasta
printf '@r1\nACGT\nAC\n+\nIIII\nII\n' > wrapped.fastq
printf '@r1\nACGT\n+\nII\n@r2\nAC\n+\nII\n' > short_quality.fastq
printf 'junk\n@r1\nACGT\n+\nIIII\n' > leading_junk.fastq

cat > check.py <<'EOF'
import importlib.util
import shutil
import sys

spec = importlib.util.spec_from_file_location('wrapper', '/repo_bin/run_microbe_census_py3fix.py')
wrapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wrapper)

passed = failed = 0
def check(desc, ok):
    global passed, failed
    print(('✓ ' if ok else '✗ ') + desc)
    passed, failed = passed + ok, failed + (not ok)

wrapper.check_version(shutil.which('run_microbe_census.py'))
check('version guard accepts the installed v1.1.1 code', True)

def same(files, threads=4):
    new = wrapper.count_bases({'verbose': False, 'seqfiles': files, 'threads': threads})
    old = wrapper.stock_count_bases({'verbose': False, 'seqfiles': files})
    return new == old, new, old

for name, fast in [('plain.fastq', True), ('plain.fastq.gz', True), ('no_last_newline.fastq', True),
                   ('empty_sequence.fastq', True), ('quality_at.fastq', True), ('empty.fastq', True),
                   ('crlf.fastq', False), ('seqs.fasta', False), ('wrapped.fastq', False),
                   ('short_quality.fastq', False), ('leading_junk.fastq', False)]:
    ok, new, old = same([name])
    path = 'fast path' if fast else 'stock fallback'
    check(f'{name}: {new} = stock {old} ({path})',
          ok and ((wrapper.fastq_bases(name) is not None) == fast))

files = ['plain.fastq.gz', 'crlf.fastq', 'seqs.fasta', 'no_last_newline.fastq']
for threads in (1, 2, 8):
    ok, new, old = same(files, threads)
    check(f'four files, {threads} worker(s): {new} = stock {old}', ok)

mc = wrapper.mc
original = mc.read_seqfile
def read_seqfile(infile):
    return original(infile)
mc.read_seqfile = read_seqfile
try:
    wrapper.check_version(shutil.which('run_microbe_census.py'))
    check('version guard refuses a changed read_seqfile', False)
except SystemExit as stop:
    check('version guard refuses a changed read_seqfile', 'read_seqfile' in str(stop))
mc.read_seqfile = original

print(f'\n{passed} passed, {failed} failed')
sys.exit(1 if failed else 0)
EOF

docker run --rm -u "$(id -u):$(id -g)" \
    -v "${TMP_DIR}:${TMP_DIR}" -w "${TMP_DIR}" \
    -v "${REPO_DIR}/bin:/repo_bin:ro" \
    "${MC_IMG}" python check.py
