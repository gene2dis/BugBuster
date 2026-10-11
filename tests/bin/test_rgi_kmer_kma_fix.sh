#!/bin/bash
#
# Regression tests for bin/rgi_kmer_kma_fix.py (RGI_KMER). RGI 6.0.3's
# CARDkmers.get_bwt_sequences read KMA BAMs from the wrong SAM columns (RNAME
# as the flag, QUAL as the sequence), so every kmer_query worker exited and
# the parent hung. Runs inside the pinned RGI image (the patch imports RGI's
# own module) and checks:
#   - the version guard accepts the installed 6.0.3 method and refuses any
#     other one;
#   - on a BAM whose @PG is KMA, the unpatched method writes ARO names as
#     flags (the bug) and the patched one writes numeric flags and the real
#     read sequences;
#   - a worker that dies without a result (out-of-memory kill, or exit()
#     without a result) makes the patched execute_threads exit with the
#     worker's status, where the unpatched one waits forever;
#   - the time limit stops the worker processes and exits with status 124.
#
# Requirements: docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${SCRIPT_DIR}/../.." && pwd )"
RGI_IMG="quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0"

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

echo "=== rgi_kmer_kma_fix.py test suite ==="

# Three reads aligned to a CARD-like reference, @PG written by KMA, as RGI_BWT
# produces with --rgi_aligner kma. RGI's extraction loop (`while read line;
# do awk ...; done`) consumes the first record before awk starts, so two of
# the three records reach the FASTA; the patch keeps that upstream behaviour.
REF='ARO:3003198|ID:41|Name:rmtH|NCBI:KC544262.1'
SEQ1="$(printf 'ACGT%.0s' {1..10})"
SEQ2="$(printf 'TTGCA%.0s' {1..8})"
SEQ3="$(printf 'GGCAT%.0s' {1..8})"
QUAL="$(printf 'I%.0s' {1..40})"
cat > "${TMP_DIR}/kma.sam" <<EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:${REF}	LN:1000
@PG	ID:KMA	PN:kma	VN:1.4.9
read_1/1	65	${REF}	1	60	40M	*	0	0	${SEQ1}	${QUAL}
read_2/1	65	${REF}	11	42	40M	*	0	0	${SEQ2}	${QUAL}
read_3/2	129	${REF}	21	30	40M	*	0	0	${SEQ3}	${QUAL}
EOF
in_rgi samtools view -b -o kma.bam kma.sam

cat > "${TMP_DIR}/run_tests.py" <<'PY'
import importlib.util
import multiprocessing
import os
import subprocess
import sys
import time

shim_path, work = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location('rgi_kmer_kma_fix', shim_path)
shim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shim)

from app.kmer_query import CARDkmers

results = []
def check(desc, ok):
    results.append(ok)
    print(('PASS ' if ok else 'FAIL ') + desc)

def extract(method, name):
    obj = object.__new__(CARDkmers)
    obj.input_bam_file = os.path.join(work, 'kma.bam')
    obj.fasta_file = os.path.join(work, name)
    method(obj)
    with open(obj.fasta_file) as fh:
        lines = fh.read().splitlines()
    return [l[1:].split('__') for l in lines[0::2]], lines[1::2]

# Unpatched 6.0.3 method on a KMA BAM: the bug
headers, seqs = extract(CARDkmers.get_bwt_sequences, 'orig.fasta')
check('unpatched method writes the reference name as the flag on a KMA BAM',
      len(headers) == 2 and all(not h[2].isdigit() for h in headers))

# Guard: refuses a changed method, then patches the installed one
class Changed(object):
    def get_bwt_sequences(self):
        return None
    def execute_threads(self):
        return None
try:
    shim.patch(Changed)
    check('guard refuses a changed method', False)
except SystemExit as e:
    check('guard refuses a changed method', 'not the 6.0.3 method' in str(e))
shim.patch(CARDkmers)
check('guard accepts the installed RGI 6.0.3 method',
      CARDkmers.get_bwt_sequences is shim.get_bwt_sequences
      and CARDkmers.execute_threads is shim.execute_threads)

headers, seqs = extract(CARDkmers.get_bwt_sequences, 'fixed.fasta')
check('patched method: qname__rname__flag__mapq headers',
      [h for h in headers] == [
          ['read_2/1', 'ARO:3003198|ID:41|Name:rmtH|NCBI:KC544262.1', '65', '42'],
          ['read_3/2', 'ARO:3003198|ID:41|Name:rmtH|NCBI:KC544262.1', '129', '30']])
check('patched method: read sequences, not qualities',
      seqs == ['TTGCA' * 8, 'GGCAT' * 8])
obj = object.__new__(CARDkmers)
obj.input_bam_file = os.path.join(work, 'kma.bam')
read, model, flag, mapq = obj.get_bwt_alignment_data('__'.join(headers[1]))
check('RGI parses the patched header (mate 2 read name)', read == 'read_3/2/2' and flag == '129')

# Dead workers: run execute_threads in a child interpreter on fake workers
# ('ok' sends its chunk back, 'kill' is SIGKILLed like an OOM kill, 'quiet'
# exits 0 without a result like the exit() on a parse error)
fanout = (
    'import importlib.util, os, signal, sys\n'
    f'spec = importlib.util.spec_from_file_location("s", {shim_path!r})\n'
    's = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)\n'
    'from app.kmer_query import CARDkmers\n'
    'def work(k, j, amr, chunk, type, out):\n'
    '    if chunk == "kill": os.kill(os.getpid(), signal.SIGKILL)\n'
    '    if chunk == "quiet": sys.exit()\n'
    '    out.put(chunk)\n'
    'obj = object.__new__(CARDkmers); obj.k = 61; obj.query_sequences = work\n'
    'method = s.execute_threads if sys.argv[1] == "patched" else CARDkmers.execute_threads\n'
    'res = method(obj, sys.argv[2].split(","), {}, None, "bwt")\n'
    'print(",".join(sorted(res)))\n'
)
def fan(which, chunks, timeout=60):
    return subprocess.run([sys.executable, '-B', '-c', fanout, which, chunks],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, timeout=timeout)
try:
    fan('original', 'ok,kill', timeout=20)
    check('unpatched execute_threads waits forever on a killed worker', False)
except subprocess.TimeoutExpired:
    check('unpatched execute_threads waits forever on a killed worker', True)
proc = fan('patched', 'a,b,c')
check('patched execute_threads returns every worker result',
      proc.returncode == 0 and proc.stdout.strip() == 'a,b,c')
proc = fan('patched', 'ok,kill,ok')
check('killed worker: exit 137 (Nextflow retries it with more memory)',
      proc.returncode == 137 and 'killed by signal 9' in proc.stderr)
proc = fan('patched', 'ok,quiet')
# (Python 3.8 gives a bare exit() in a worker exit code 1, not 0)
check('worker without a result: exit 1 with a message',
      proc.returncode == 1 and ('before sending its result' in proc.stderr
                                or 'without a result' in proc.stderr))

# Time limit: a stalled run is stopped with status 124 and no workers left
probe = (
    'import importlib.util, multiprocessing, time, sys\n'
    f'spec = importlib.util.spec_from_file_location("s", {shim_path!r})\n'
    's = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)\n'
    'p = multiprocessing.Process(target=time.sleep, args=(120,)); p.start()\n'
    'open(sys.argv[1], "w").write(str(p.pid))\n'
    's.arm_time_limit(2)\n'
    'p.join()\n'
)
pid_file = os.path.join(work, 'worker.pid')
start = time.time()
proc = subprocess.run([sys.executable, '-B', '-c', probe, pid_file],
                      stderr=subprocess.PIPE, universal_newlines=True)
elapsed = time.time() - start
check('time limit exits with status 124', proc.returncode == 124)
check('time limit fires after ~2 s', 1.5 < elapsed < 20)
check('time limit message names the override', 'ext.time_limit' in proc.stderr)
worker = int(open(pid_file).read())
time.sleep(0.5)
try:
    os.kill(worker, 0)
    alive = not open(f'/proc/{worker}/stat').read().split(') ')[1].startswith('Z')
except (ProcessLookupError, FileNotFoundError):
    alive = False
check('time limit stops the worker process', not alive)

# CLI: rejects a missing --time-limit
proc = subprocess.run([sys.executable, '-B', shim_path, 'kmer_query', '--help'],
                      stderr=subprocess.PIPE, universal_newlines=True)
check('CLI requires --time-limit', proc.returncode != 0 and 'usage' in proc.stderr)

sys.exit(0 if all(results) else 1)
PY

while IFS= read -r line; do
    case "${line}" in
        PASS\ *) check "${line#PASS }" true ;;
        FAIL\ *) check "${line#FAIL }" false ;;
    esac
done < <(in_rgi python3 -B run_tests.py "${REPO_DIR}/bin/rgi_kmer_kma_fix.py" "${TMP_DIR}" 2>/dev/null || true)

check "all 15 in-container checks reported" test $((PASS + FAIL)) -eq 15

echo ""
echo "=== ${PASS} passed, ${FAIL} failed ==="
[ "${FAIL}" -eq 0 ]
