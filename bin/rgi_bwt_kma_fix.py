#!/usr/bin/env python3
"""Delegating wrapper around RGI's `rgi` CLI for RGI_BWT (rgi bwt --aligner kma).

Usage: rgi_bwt_kma_fix.py <rgi arguments...>     e.g. rgi_bwt_kma_fix.py bwt --read_one ...

RGI 6.0.3 BWT.align_kma / align_kma_interleaved write KMA's SAM output for
every read, although well under 1% of a metagenome maps to CARD/WildCARD.
Every later step of rgi bwt works on that file: two BAM copies, a read
table, and then, for each reference sequence with reads (thousands with
WildCARD), one scan of the whole read table. On an 81M-read sample that is
81M rows scanned 7,859 times, hours per task, and the files hold ~390 bytes
per clean read on disk.

This wrapper replaces both methods with copies that run the same KMA command
and keep, of its SAM output, the header and the records placed on a
reference (RNAME not '*' and POS not 0). The others are unplaced in the BAM
rgi builds (KMA writes reads with k-mer hits but no alignment at POS 0, which
samtools turns into unplaced unmapped reads), so they reach no rgi bwt
output used downstream: the BAM's placed reads, and so the gene and allele
tables, are the same. Only overall_mapping_stats.txt and
artifacts_mapping_stats.txt change: their read totals now count the placed
reads only. The copies also stop when KMA fails; the originals ignore it, so
a failed alignment looked like a sample without hits.

KMA's own outputs (consensus, VCF: the SNPs column of the allele table) vary
between identical KMA runs with several threads, with or without this
wrapper.

The patches apply only to the exact RGI 6.0.3 methods: if an installed
source differs (an upstream fix or any other change), the wrapper stops with
an error instead of patching blindly. Update or delete it when RGI is
re-pinned. Other aligners (--aligner bowtie2|bwa) run unchanged.
"""

import hashlib
import inspect
import os
import runpy
import shutil
import subprocess
import sys

# sha256 of inspect.getsource(BWT.<method>) as shipped in RGI 6.0.3
# (quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0)
RGI_603_METHOD_SHA256 = {
    'align_kma':             '929c984dbb6641ef11722365c07115cfeba9ffa1f80b9225f9056b9befa14e2a',
    'align_kma_interleaved': '982060f17da7819cac1cab83e33d1004e13de3f89e211e83a04ae9a0c2f6d10b',
}

# SAM header lines and the records placed on a reference
PLACED_READS_AWK = '/^@/ || ($3 != "*" && $4 != 0)'

# KMA 1.4.9 (the pinned image's) can exit 95 after a successful run; real
# errors give other statuses (2 for a missing index, 137 for a kill)
KMA_OK_STATUS = (0, 95)

logger = None   # RGI's own logger (app.BWT.logger), set in main()


def run_kma(cmd, output_sam_file):
    """Run a KMA command line (SAM on stdout) and keep the placed reads.

    Success needs an accepted KMA exit status, a successful filter, and
    KMA's result table (<output_sam_file>.temp.res, written when KMA closes
    its files; rgi bwt reads it later for the read coverage).
    """
    with open(output_sam_file, 'w') as out:
        kma = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        awk = subprocess.Popen(['awk', '-F', '\t', PLACED_READS_AWK], stdin=kma.stdout, stdout=out)
        kma.stdout.close()   # awk holds the pipe now
        awk_status, kma_status = awk.wait(), kma.wait()
    if kma_status not in KMA_OK_STATUS or not os.path.exists(f'{output_sam_file}.temp.res'):
        sys.exit(f'rgi_bwt_kma_fix.py: KMA failed (exit status {kma_status}): {cmd}')
    if awk_status != 0:
        sys.exit(f'rgi_bwt_kma_fix.py: filtering the KMA SAM failed (awk exit status {awk_status})')


def align_kma(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_kma (RGI 6.0.3): placed reads only, KMA exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    logger.info("align reads -ipe {} {} to {}".format(self.read_one, self.read_two, reference_genome))
    run_kma("kma -mem_mode -ex_mode -1t1 -vcf -ipe {read_one} {read_two} -t {threads} -t_db {index_directory} "
            "-o {output_sam_file}.temp -sam".format(
                threads=self.threads, index_directory=index_directory, read_one=self.read_one,
                read_two=self.read_two, output_sam_file=output_sam_file), output_sam_file)


def align_kma_interleaved(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_kma_interleaved (RGI 6.0.3): placed reads only, KMA exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    logger.info("align reads -int {} to {}".format(self.read_one, reference_genome))
    run_kma("kma -mem_mode -ex_mode -1t1 -vcf -int {read_one} -t {threads} -t_db {index_directory} "
            "-o {output_sam_file}.temp -sam".format(
                threads=self.threads, index_directory=index_directory, read_one=self.read_one,
                output_sam_file=output_sam_file), output_sam_file)


def patch(bwt_cls):
    """Replace the BWT methods after checking they are the 6.0.3 ones."""
    replacements = (('align_kma', align_kma), ('align_kma_interleaved', align_kma_interleaved))
    for name, _ in replacements:
        try:
            installed = inspect.getsource(getattr(bwt_cls, name))
        except (OSError, TypeError, AttributeError):
            installed = ''   # no source to verify: treat as not the 6.0.3 method
        digest = hashlib.sha256(installed.encode()).hexdigest()
        if digest != RGI_603_METHOD_SHA256[name]:
            sys.exit(
                f'rgi_bwt_kma_fix.py: the installed RGI BWT.{name} is not the '
                f'6.0.3 method this patch targets (sha256 {digest}). Check '
                'whether this RGI version still writes every read to the KMA '
                'SAM, then update or delete bin/rgi_bwt_kma_fix.py.')
    for name, replacement in replacements:
        setattr(bwt_cls, name, replacement)


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit('usage: rgi_bwt_kma_fix.py <rgi arguments...>')

    global logger
    import app.BWT
    logger = app.BWT.logger
    patch(app.BWT.BWT)

    rgi = shutil.which('rgi')
    if rgi is None:
        sys.exit('rgi_bwt_kma_fix.py: rgi not found on PATH')
    sys.argv = [rgi] + args
    runpy.run_path(rgi, run_name='__main__')


if __name__ == '__main__':
    main()
