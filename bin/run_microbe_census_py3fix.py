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

This wrapper replaces check_rapsearch with the corrected implementation
(a faithful copy of the upstream fix, tolerant of str or bytes) and then
executes the stock CLI unchanged, forwarding all arguments. Delete it if
a fixed MicrobeCensus build ever reaches bioconda.
"""

import runpy
import shutil
import sys

from microbe_census import microbe_census as mc


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


def main():
    mc.check_rapsearch = check_rapsearch

    cli = shutil.which('run_microbe_census.py')
    if cli is None:
        sys.exit('run_microbe_census.py not found on PATH — is this running '
                 'inside the microbecensus container?')
    sys.argv[0] = cli
    runpy.run_path(cli, run_name='__main__')


if __name__ == '__main__':
    main()
