#!/usr/bin/env python3
"""Delegating wrapper around RGI's `rgi` CLI for RGI_KMER (rgi kmer_query --bwt).

Usage: rgi_kmer_kma_fix.py --time-limit SECONDS <rgi arguments...>

1. KMA BAMs. RGI 6.0.3 CARDkmers.get_bwt_sequences reads the aligner from the
   BAM's last @PG line and, for KMA ("ID:KMA"), extracts the alignment fields
   from the wrong SAM columns: the reference name (column 3) as the flag and
   the base qualities (column 11) as the read sequence. KMA (1.4.9 in the
   pinned image) writes standard SAM. Every kmer_query worker then logs
   "failed to parse BAM file" and exits, and the parent waits forever: with
   KMA, RGI's default aligner and ours, RGI_KMER hung until its time limit.
   This wrapper replaces the method with a copy that uses the standard
   column order (the one RGI already uses for bowtie2/bwa BAMs) for every
   aligner; the extraction command is otherwise unchanged.

2. Dead workers. CARDkmers.execute_threads forks one worker per thread and
   then blocks on one queue result per worker. A worker that dies without
   sending one (killed for exceeding the task memory, or the exit() above)
   leaves the parent waiting forever, so the task never fails and Nextflow
   never retries it with more memory. The replacement waits the same way
   but checks the workers: if one has died without a result, it stops the
   others and exits with that worker's status (128 + signal for a killed
   worker, e.g. 137 for an out-of-memory kill).

3. Time limit. If rgi has not finished after SECONDS, the wrapper stops its
   worker processes and exits with status 124 and a clear error, so any
   other stall fails fast instead of holding the task until Nextflow's own
   time limit (which, on the local executor, aborts the whole run).

The patches apply only to the exact RGI 6.0.3 methods: if an installed
source differs (an upstream fix or any other change), the wrapper stops with
an error instead of patching blindly. Update or delete it when RGI is
re-pinned.
"""

import hashlib
import inspect
import multiprocessing
import os
import queue
import runpy
import shutil
import signal
import sys

# sha256 of inspect.getsource(CARDkmers.<method>) as shipped in RGI 6.0.3
# (quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0)
RGI_603_METHOD_SHA256 = {
    'get_bwt_sequences': '0316d8dbd9faa426747898be43fb5ffd6d609b593fd651f825a3cb2b98f48bd0',
    'execute_threads':   '854cd389be28060eeb5164beb8ca82824bc4da7af5b528d9640a1681db64c045',
}

TIMEOUT_EXIT_STATUS = 124
WORKER_POLL_SECONDS = 5


def get_bwt_sequences(self):
    """Fixed copy of CARDkmers.get_bwt_sequences (RGI 6.0.3).

    Filter out unmapped and supplementary alignments; store QNAME, RNAME,
    FLAG and MAPQ in the FASTA header and SEQ as the sequence. The 6.0.3
    method used columns 1/4/3/5 and 11 for KMA BAMs; every aligner now
    takes the standard columns 1/3/2/5 and 10.
    """
    os.system("samtools index {input_bam}".format(input_bam=self.input_bam_file))
    os.system("""samtools view -F 4 -F 2048 {bam} | while read line; do awk '{cmd}'; done > {out}"""
                .format(bam=self.input_bam_file, cmd="""{print ">"$1"__"$3"__"$2"__"$5"\\n"$10}""", out=self.fasta_file))


def stop_workers(processes, status, message):
    """Kill the still-running workers, report, and exit with `status`."""
    for process in processes:
        if process.is_alive():
            process.kill()
    sys.stderr.write(f'rgi_kmer_kma_fix.py: {message}\n')
    sys.stderr.flush()
    os._exit(status)


def execute_threads(self, split_sequences, j, amr_kmers, type):
    """Fixed copy of CARDkmers.execute_threads (RGI 6.0.3).

    Same fan-out and result order as the original; instead of blocking on
    the queue forever, it polls and fails as soon as a worker has died
    without sending its result.
    """
    output = multiprocessing.Queue()
    processes = []
    for ind in range(len(split_sequences)):
        process = multiprocessing.Process(target=self.query_sequences,
            args=(self.k, j, amr_kmers, split_sequences[ind], type, output))
        process.start()
        processes.append(process)
    results = []
    while len(results) < len(processes):
        try:
            results.append(output.get(timeout=WORKER_POLL_SECONDS))
            continue
        except queue.Empty:
            pass
        failed = [p for p in processes if p.exitcode not in (None, 0)]
        if failed:
            code = failed[0].exitcode
            status = 128 - code if code < 0 else code
            reason = (f'was killed by signal {-code} (out of memory?)'
                      if code < 0 else f'exited with status {code}')
            stop_workers(processes, status,
                         f'an rgi kmer_query worker {reason} before sending its '
                         f'result; stopped the run (exit {status}).')
        if all(p.exitcode is not None for p in processes):
            # Every worker has exited: anything they sent is in the pipe now
            try:
                while len(results) < len(processes):
                    results.append(output.get(timeout=WORKER_POLL_SECONDS))
            except queue.Empty:
                stop_workers(processes, 1,
                             f'{len(processes) - len(results)} rgi kmer_query '
                             'worker(s) exited without a result (see the rgi '
                             'errors above); stopped the run.')
    for process in processes:
        process.join()
    return results


def patch(kmers_cls):
    """Replace the CARDkmers methods after checking they are the 6.0.3 ones."""
    for name, replacement in (('get_bwt_sequences', get_bwt_sequences),
                              ('execute_threads', execute_threads)):
        try:
            installed = inspect.getsource(getattr(kmers_cls, name))
        except (OSError, TypeError, AttributeError):
            installed = ''   # no source to verify: treat as not the 6.0.3 method
        digest = hashlib.sha256(installed.encode()).hexdigest()
        if digest != RGI_603_METHOD_SHA256[name]:
            sys.exit(
                f'rgi_kmer_kma_fix.py: the installed RGI CARDkmers.{name} is '
                f'not the 6.0.3 method this patch targets (sha256 {digest}). '
                'Check whether this RGI version still needs the fix, then '
                'update or delete bin/rgi_kmer_kma_fix.py.')
    for name, replacement in (('get_bwt_sequences', get_bwt_sequences),
                              ('execute_threads', execute_threads)):
        setattr(kmers_cls, name, replacement)


def arm_time_limit(seconds):
    """After `seconds`, stop all child processes and exit with status 124."""
    def on_timeout(signum, frame):
        for child in multiprocessing.active_children():
            child.kill()
        sys.stderr.write(
            f'rgi_kmer_kma_fix.py: rgi did not finish within {seconds} s; '
            'stopped it. Raise the limit with ext.time_limit for RGI_KMER '
            '(docs/troubleshooting.md) if this sample is just large.\n')
        sys.stderr.flush()
        os._exit(TIMEOUT_EXIT_STATUS)
    signal.signal(signal.SIGALRM, on_timeout)
    signal.alarm(seconds)


def main():
    args = sys.argv[1:]
    if len(args) < 2 or args[0] != '--time-limit' or not args[1].isdigit():
        sys.exit('usage: rgi_kmer_kma_fix.py --time-limit SECONDS <rgi arguments...>')
    seconds = int(args[1])

    from app.kmer_query import CARDkmers
    patch(CARDkmers)

    rgi = shutil.which('rgi')
    if rgi is None:
        sys.exit('rgi_kmer_kma_fix.py: rgi not found on PATH')
    if seconds > 0:
        arm_time_limit(seconds)
    sys.argv = [rgi] + args[2:]
    runpy.run_path(rgi, run_name='__main__')


if __name__ == '__main__':
    main()
