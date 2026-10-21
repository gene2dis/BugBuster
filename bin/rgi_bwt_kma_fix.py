#!/usr/bin/env python3
"""Delegating wrapper around RGI's `rgi` CLI for RGI_BWT (rgi bwt --aligner kma).

Usage: rgi_bwt_kma_fix.py <rgi arguments...>     e.g. rgi_bwt_kma_fix.py bwt --read_one ...

RGI 6.0.8 BWT.align_kma / align_kma_interleaved write KMA's SAM output for
every read, although well under 1% of a metagenome maps to CARD/WildCARD.
Every later step of rgi bwt works on that file: samtools fixmate, two BAM
copies, and a read table that rgi then loads into memory (6.0.3 scanned the
whole table once per reference with reads instead, thousands of passes with
WildCARD). Unfiltered, the files hold ~390 bytes per clean read on disk and
the table would hold every read in memory.

This wrapper replaces both methods with copies that run the same KMA command
and keep, of its SAM output, the header and the records placed on a
reference (RNAME not '*' and POS not 0), grouped by read name. KMA writes
reads with k-mer hits but no alignment at POS 0 and reads without hits with
RNAME '*'; neither is placed in the BAM. The copies also stop when KMA
fails; the originals ignore it, so a failed alignment looked like a sample
without hits.

Grouping by name makes the BAM deterministic. 6.0.8 runs `samtools fixmate
-z all -m` on the SAM, which pairs only mates on adjacent lines, and KMA
writes its records in a thread-dependent order (mates are often apart, and
the order changes between identical runs). On adjacent mates, fixmate sets
the mate fields and moves an unmapped mate to its partner's position, where
idxstats counts it as a read "with flanking sequence" (part of "All Mapped
Reads"). Unfiltered, those counts and the allele table's Mate Pair Linkage
change from run to run. With the filter, every pair whose two mates are
placed gets its mate fields, and an unmapped mate is never counted, as in
6.0.3 (flanking reads 0). overall_mapping_stats.txt and
artifacts_mapping_stats.txt count the placed reads only.

KMA's own outputs (consensus, VCF: the SNPs column of the allele table) vary
between identical KMA runs with several threads, with or without this
wrapper.

The patches apply only to the exact RGI 6.0.8 methods: if an installed
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

# sha256 of inspect.getsource(BWT.<method>) as shipped in RGI 6.0.8
# (quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0)
RGI_608_METHOD_SHA256 = {
    'align_kma':             '6ff96153fa4bdfbd2f563e163c99004593a33b2821064b82255336c1e868e7ec',
    'align_kma_interleaved': '6406e3246ebc6c4ba1482b595bdeee4046c139ce40da3506d7a6c15bb708a745',
}

# SAM header lines and the records placed on a reference
PLACED_READS_AWK = '/^@/ || ($3 != "*" && $4 != 0)'

# KMA 1.4.9 (RGI 6.0.3's image) could exit 95 after a successful run; the
# 1.6.11 of the pinned image exited 0 in every test, and 95 is still accepted
# when the result table exists. Real errors give other statuses (2 for a
# missing index, 137 for a kill)
KMA_OK_STATUS = (0, 95)

logger = None   # RGI's own logger (app.BWT.logger), set in main()


def run_kma(cmd, output_sam_file):
    """Run a KMA command line (SAM on stdout) and keep the placed reads.

    Success needs an accepted KMA exit status, a successful filter, and
    KMA's result table (<output_sam_file>.temp.res, written when KMA closes
    its files; rgi bwt reads it later for the read coverage). The placed
    records (well under 1% of the reads) are then grouped by read name.
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
    group_by_name(output_sam_file)


def group_by_name(sam_file):
    """Rewrite a SAM file with its records stably sorted by QNAME (header first).

    The sort runs in child processes, so its memory (about 1 GB for the
    3 M placed records of an ARG-rich deep sample) is released before rgi
    loads its databases.
    """
    grouped = f'{sam_file}.byname'
    with open(grouped, 'wb') as out:
        header = subprocess.run(['awk', '/^@/', sam_file], stdout=out)
        out.flush()
        records = subprocess.Popen(['awk', '!/^@/', sam_file], stdout=subprocess.PIPE)
        by_name = subprocess.Popen(['sort', '-s', '-t', '\t', '-k1,1'], stdin=records.stdout, stdout=out,
                                   env=dict(os.environ, LC_ALL='C'))
        records.stdout.close()   # sort holds the pipe now
        sort_status, records_status = by_name.wait(), records.wait()
    if header.returncode or records_status or sort_status:
        sys.exit(f'rgi_bwt_kma_fix.py: grouping the KMA SAM by read name failed (awk {header.returncode}, '
                 f'{records_status}; sort {sort_status})')
    os.replace(grouped, sam_file)


def align_kma(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_kma (RGI 6.0.8): placed reads by name, KMA exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    logger.info("align reads -ipe {} {} to {}".format(self.read_one, self.read_two, reference_genome))
    run_kma("kma -mem_mode -ex_mode -1t1 -vcf -ipe {read_one} {read_two} -t {threads} -t_db {index_directory} "
            "-o {output_sam_file}.temp -sam".format(
                threads=self.threads, index_directory=index_directory, read_one=self.read_one,
                read_two=self.read_two, output_sam_file=output_sam_file), output_sam_file)


def align_kma_interleaved(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_kma_interleaved (RGI 6.0.8): placed reads by name, KMA exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    logger.info("align reads -int {} to {}".format(self.read_one, reference_genome))
    run_kma("kma -mem_mode -ex_mode -1t1 -vcf -int {read_one} -t {threads} -t_db {index_directory} "
            "-o {output_sam_file}.temp -sam".format(
                threads=self.threads, index_directory=index_directory, read_one=self.read_one,
                output_sam_file=output_sam_file), output_sam_file)


def patch(bwt_cls):
    """Replace the BWT methods after checking they are the 6.0.8 ones."""
    replacements = (('align_kma', align_kma), ('align_kma_interleaved', align_kma_interleaved))
    for name, _ in replacements:
        try:
            installed = inspect.getsource(getattr(bwt_cls, name))
        except (OSError, TypeError, AttributeError):
            installed = ''   # no source to verify: treat as not the 6.0.8 method
        digest = hashlib.sha256(installed.encode()).hexdigest()
        if digest != RGI_608_METHOD_SHA256[name]:
            sys.exit(
                f'rgi_bwt_kma_fix.py: the installed RGI BWT.{name} is not the '
                f'6.0.8 method this patch targets (sha256 {digest}). Check '
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
