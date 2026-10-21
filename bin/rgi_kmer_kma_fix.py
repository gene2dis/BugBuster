#!/usr/bin/env python3
"""Delegating wrapper around RGI's `rgi` CLI for RGI_KMER (rgi kmer_query --bwt).

Usage: rgi_kmer_kma_fix.py --time-limit SECONDS <rgi arguments...>

1. Dead workers. CARDkmers.execute_threads forks one worker per thread and
   then blocks on one queue result per worker. A worker that dies without
   sending one (killed for exceeding the task memory, or exiting on an
   error) leaves the parent waiting forever, so the task never fails and
   Nextflow never retries it with more memory. The replacement waits the same way
   but checks the workers: if one has died without a result, it stops the
   others and exits with that worker's status (128 + signal for a killed
   worker, e.g. 137 for an out-of-memory kill).

2. Time limit. If rgi has not finished after SECONDS, the wrapper stops its
   worker processes and exits with status 124 and a clear error, so any
   other stall fails fast instead of holding the task until Nextflow's own
   time limit (which, on the local executor, aborts the whole run).

RGI 6.0.3 also read KMA BAMs from the wrong SAM columns, so every worker
failed; this wrapper patched that too until RGI 6.0.8, whose
CARDkmers.get_bwt_sequences reads the BAM with pysam (checked on real BAMs:
the same FASTA, plus the first alignment, which the 6.0.3 shell loop
dropped).

The patch applies only to the exact RGI 6.0.8 method: if the installed
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

# sha256 of inspect.getsource(CARDkmers.execute_threads) as shipped in RGI
# 6.0.8 (quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0)
RGI_608_METHOD_SHA256 = {
    'execute_threads': '8a8664bf2a1205be5141cc872588dc088ed2a0f4ed02482901bdc8f9ec699be4',
}

TIMEOUT_EXIT_STATUS = 124
WORKER_POLL_SECONDS = 5


def stop_workers(processes, status, message):
    """Kill the still-running workers, report, and exit with `status`."""
    for process in processes:
        if process.is_alive():
            process.kill()
    sys.stderr.write(f'rgi_kmer_kma_fix.py: {message}\n')
    sys.stderr.flush()
    os._exit(status)


def execute_threads(self, split_sequences, j, amr_kmers, type):
    """Fixed copy of CARDkmers.execute_threads (RGI 6.0.8).

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
    """Replace the CARDkmers methods after checking they are the 6.0.8 ones."""
    for name in RGI_608_METHOD_SHA256:
        try:
            installed = inspect.getsource(getattr(kmers_cls, name))
        except (OSError, TypeError, AttributeError):
            installed = ''   # no source to verify: treat as not the 6.0.8 method
        digest = hashlib.sha256(installed.encode()).hexdigest()
        if digest != RGI_608_METHOD_SHA256[name]:
            sys.exit(
                f'rgi_kmer_kma_fix.py: the installed RGI CARDkmers.{name} is '
                f'not the 6.0.8 method this patch targets (sha256 {digest}). '
                'Check whether this RGI version still needs the fix, then '
                'update or delete bin/rgi_kmer_kma_fix.py.')
    kmers_cls.execute_threads = execute_threads


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
