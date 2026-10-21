#!/usr/bin/env python3
"""Delegating wrapper around MicrobeCensus' run_microbe_census.py (design doc Q4, T5).

The pinned biocontainer (microbecensus:1.1.1--pyhca03a8a_2, the only
Python-3 build) ships the v1.1.1 tarball, which predates the upstream
Python-3 fixes: check_rapsearch() calls .split('\\n') on the bytes stderr
of Popen.communicate() and crashes before any reads are processed
(upstream issue #32; the fix, commit 32eb3620, was never released — the
project is unmaintained). The Python-2 biocontainer is broken differently
(its Buildroot base lacks the libstdc++ needed by the bundled RAPsearch2
binary), so shimming the py3 image is the only working configuration.

This wrapper replaces two functions and then executes the stock CLI
unchanged, forwarding all arguments:

- check_rapsearch: the corrected implementation (a faithful copy of the
  upstream fix, tolerant of str or bytes);
- count_bases: the genome-equivalents step counts every base of every
  input file in one Python loop, about 16 Mbp/s on one core, which at full
  depth (~18 Gbp) took about as long as the average-genome-size estimate
  (performance review M6). The replacement counts the files in parallel
  (one process per file, up to -t) and reads plain 4-line FASTQ in large
  blocks; any file that is not plain 4-line FASTQ (FASTA, wrapped or
  malformed records, CR line ends, non-ASCII, bzip2) is counted by the
  stock function, so total_bases is the stock value in every case.

The patches apply only to the exact v1.1.1 code: if the installed source
of a patched function, of a function the stock count depends on, or of the
CLI differs, the wrapper stops instead of patching. Delete it if a fixed
MicrobeCensus build ever reaches bioconda.
"""

import gzip
import hashlib
import inspect
import operator
import runpy
import shutil
import sys
from itertools import repeat
import multiprocessing

from microbe_census import microbe_census as mc

# sha256 of the installed v1.1.1 source (inspect.getsource for functions,
# the file for the CLI) in microbecensus:1.1.1--pyhca03a8a_2
V111_SHA256 = {
    'check_rapsearch': 'b7aa269483cc01d3ba1086084678dfd6ecf28d0222a3c2a91e2ed214b6e2b18c',
    'count_bases': '7960d0e88054abe10460de62339692eae84c698f7764db42f47f386cb4cb5544',
    'read_seqfile': 'd5877c51c219be62bc88d5796cbab8adf8933360ff375dc521dfc776103173b0',
    'open_file': '9f72f1386e8c78cb1ea7bec4cd36c3245b2504485ffcf0445c6bffaea3975c3b',
    'cli': '3a98134604c4a848f732639379c0913f34a5df182bc5704c86f08dd83366d688',
}

CHUNK = 8 * 1024 * 1024
SEQ_LINE_STARTS = (b'@', b'+', b'>')
stock_count_bases = mc.count_bases


def check_rapsearch(rapsearch):
    """Corrected copy of microbe_census.check_rapsearch (decodes stderr)."""
    import os
    import subprocess

    process = subprocess.Popen(
        rapsearch + ' -h', shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    _output, error = process.communicate()
    if isinstance(error, bytes):
        error = error.decode()
    if os.path.isdir(rapsearch):
        sys.exit("Problem executing rapsearch2: '%s'" % rapsearch)
    elif len(error.split('\n')) < 2:
        sys.exit("Problem executing rapsearch2: '%s'" % rapsearch)
    elif error.split('\n')[1] != 'rapsearch v2.15: Fast protein similarity search tool for short reads':
        sys.exit("Incorrect version of rapsearch2 detected:'%s\nMicrobeCensus requires rapsearch v2.15" % rapsearch)


def fastq_bases(path):
    """Bases of a plain 4-line FASTQ, as the stock count gives them, or None.

    The stock parser (read_seqfile on a text-mode open_file) counts each
    record's sequence line without its newline. Records are checked so that
    the stock parser would read the same lines: '@' header, sequence line not
    starting with '@', '+' or '>', '+' line, quality at least as long as the
    sequence; no CR (text mode turns CR into line ends) and ASCII only.
    """
    ext = path.split('.')[-1]
    if ext == 'gz':
        handle = gzip.open(path, 'rb')
    elif ext == 'bz2':
        return None
    else:
        handle = open(path, 'rb')
    total, rest = 0, b''
    with handle:
        while True:
            chunk = handle.read(CHUNK)
            if chunk:
                if b'\r' in chunk or not chunk.isascii():
                    return None
                lines = (rest + chunk).split(b'\n')
                n = (len(lines) - 1) // 4 * 4
                rest = b'\n'.join(lines[n:])
            else:
                lines = rest.split(b'\n') if rest else []
                n = len(lines)
                if n % 4:
                    return None
            heads, seqs = lines[0:n:4], lines[1:n:4]
            plus, quals = lines[2:n:4], lines[3:n:4]
            if not (all(map(bytes.startswith, heads, repeat(b'@')))
                    and all(map(bytes.startswith, plus, repeat(b'+')))
                    and not any(map(bytes.startswith, seqs, repeat(SEQ_LINE_STARTS)))
                    and all(map(operator.ge, map(len, quals), map(len, seqs)))):
                return None
            total += sum(map(len, seqs))
            if not chunk:
                return total


def file_bases(index, path, results):
    try:
        bases = fastq_bases(path)
        if bases is None:
            bases = stock_count_bases({'verbose': False, 'seqfiles': [path]})
        results.put((index, bases, None))
    except BaseException as error:
        results.put((index, None, '%s: %r' % (path, error)))


def count_bases(args):
    """Parallel replacement for microbe_census.count_bases (same total).

    Forked processes (the CLI runs as __main__ under runpy, so a Pool could
    not pickle the worker function), at most -t of them at a time.
    """
    if args['verbose']:
        print("Computing number of genome equivalents...")
    ctx = multiprocessing.get_context('fork')
    files = args['seqfiles']
    width = max(1, min(len(files), args['threads']))
    results = ctx.SimpleQueue()
    total = 0
    for start in range(0, len(files), width):
        workers = [ctx.Process(target=file_bases, args=(i, path, results))
                   for i, path in enumerate(files[start:start + width], start)]
        for worker in workers:
            worker.start()
        for _ in workers:
            _index, bases, error = results.get()
            if error is not None:
                sys.exit('count_bases failed on %s' % error)
            total += bases
        for worker in workers:
            worker.join()
    return total


def check_version(cli):
    sources = {name: inspect.getsource(getattr(mc, name))
               for name in ('check_rapsearch', 'count_bases', 'read_seqfile', 'open_file')}
    with open(cli) as handle:
        sources['cli'] = handle.read()
    changed = [name for name, text in sources.items()
               if hashlib.sha256(text.encode()).hexdigest() != V111_SHA256[name]]
    if changed:
        sys.exit('run_microbe_census_py3fix.py: the installed MicrobeCensus is '
                 'not the v1.1.1 code this wrapper patches (changed: %s). '
                 'Check whether upstream fixed check_rapsearch (#32), then '
                 'update or delete bin/run_microbe_census_py3fix.py.'
                 % ', '.join(changed))


def main():
    cli = shutil.which('run_microbe_census.py')
    if cli is None:
        sys.exit('run_microbe_census.py not found on PATH — is this running '
                 'inside the microbecensus container?')
    check_version(cli)
    mc.check_rapsearch = check_rapsearch
    mc.count_bases = count_bases
    sys.argv[0] = cli
    runpy.run_path(cli, run_name='__main__')


if __name__ == '__main__':
    main()
