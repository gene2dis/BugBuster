#!/usr/bin/env python3
"""Delegating wrapper around RGI's `rgi` CLI for RGI_BWT (rgi bwt).

Usage: rgi_bwt_kma_fix.py <rgi arguments...>     e.g. rgi_bwt_kma_fix.py bwt --read_one ...

RGI 6.0.8's BWT.align_* methods write the aligner's SAM output for every
read, although well under 1% of a metagenome maps to CARD/WildCARD. Every
later step of rgi bwt works on that file: samtools fixmate, two BAM copies,
and a read table that rgi then loads into memory (6.0.3 scanned the whole
table once per reference with reads instead, thousands of passes with
WildCARD). Unfiltered, the files hold ~390 bytes per clean read on disk and
the table holds every read in memory (~1 KB per read).

This wrapper replaces those methods with copies that run the same aligner
command and filter its SAM output. The copies also stop when the aligner
fails; the originals ignore it, so a failed alignment looked like a sample
without hits.

KMA (--aligner kma, the default): the header and the records placed on a
reference (RNAME not '*' and POS not 0), grouped by read name. KMA writes
reads with k-mer hits but no alignment at POS 0 and reads without hits with
RNAME '*'; neither is placed in the BAM. Grouping by name makes the BAM
deterministic: 6.0.8 runs `samtools fixmate -z all -m` on the SAM, which
pairs only mates on adjacent lines, and KMA writes its records in a
thread-dependent order (mates are often apart, and the order changes between
identical runs). On adjacent mates, fixmate sets the mate fields and moves an
unmapped mate to its partner's position, where idxstats counts it as a read
"with flanking sequence" (part of "All Mapped Reads"). Unfiltered, those
counts and the allele table's Mate Pair Linkage change from run to run. With
the filter, every pair whose two mates are placed gets its mate fields (mates
named .../1 and .../2 stay unpaired, as fixmate needs equal names), and an
unmapped mate is never counted, as in 6.0.3 (flanking reads 0). The grouping
sort runs in a child process and holds the placed records in memory (about
1 GB per 3 M records).

bowtie2 and bwa (--aligner bowtie2|bwa): the header and every record with a
reference name, so only pairs with no mapped mate are dropped. Those reach no
rgi bwt output (they are unplaced in the BAM and counted under '*', which rgi
skips), and mates stay on adjacent lines, so the tables are those of the
unfiltered run.

overall_mapping_stats.txt and artifacts_mapping_stats.txt count the kept
records only.

A failed aligner or filter stops rgi with its exit status (128 + signal for
a killed one), so RGI_BWT retries a task killed for memory (137).

The task stages localDB/ as one symlink per file. When an index there is
incomplete, RGI rebuilds it with the aligner, which would write through
those symlinks into the shared database; the wrapper first removes the
task's symlinks to that index, so the rebuild writes new files in the task.

KMA's own outputs (consensus, VCF: the SNPs column of the allele table) vary
between identical KMA runs with several threads, with or without this
wrapper.

The patches apply only to the exact RGI 6.0.8 methods: if an installed
source differs (an upstream fix or any other change), the wrapper stops with
an error instead of patching blindly. Update or delete it when RGI is
re-pinned.
"""

import glob
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
    'align_kma':                    '6ff96153fa4bdfbd2f563e163c99004593a33b2821064b82255336c1e868e7ec',
    'align_kma_interleaved':        '6406e3246ebc6c4ba1482b595bdeee4046c139ce40da3506d7a6c15bb708a745',
    'align_bowtie2':                'abf4f424b2900d158bbac175f311cbaa3118e8d068517410aa311ddb96ac90b5',
    'align_bowtie2_unpaired':       'a4e6c2e019607c984855092b4d0dddd15f6762903913e5b441b1191eeb4ac84d',
    'align_bwa_paired_end_mapping': '4e0cee70bb8ab090f8e7d3a255c538ce335cc356e67b96c1e9115ac64d5f3d38',
    'align_bwa_single_end_mapping': 'badfde083ed212dcd63e5b3de3c28f0a66b0a238d1fcc231dd0fd8442014ec03',
    'create_index':                 'b3aafd80ea560b53c24d8200eb3bfad199f9e3de070ef22600228b772a427c3d',
}

# KMA: the header (its @HD line now says the records are grouped by query
# name) and the records placed on a reference, grouped by read name: a
# stable sort in a child process, whose output follows the flushed header.
# busybox awk's close() returns the sort's wait status (the signal number
# for a killed sort), passed on as 128 + signal
KMA_PLACED_READS_AWK = (
    'BEGIN { by_name = "LC_ALL=C sort -s -t \'\\t\' -k1,1" } '
    '/^@/ { if ($1 == "@HD") sub(/\\tGO:[a-z]+/, "\\tGO:query"); print; next } '
    '$3 != "*" && $4 != 0 { print | by_name; n++ } '
    'END { fflush(); if (n && (r = close(by_name)) != 0) exit (r % 256 ? 128 + r % 128 : 1) }'
)

# bowtie2/bwa: the header and every record with a reference name
ALIGNED_READS_AWK = '/^@/ || $3 != "*"'

# KMA 1.4.9 (RGI 6.0.3's image) could exit 95 after a successful run; the
# 1.6.11 of the pinned image exited 0 in every test, and 95 is still accepted
# when the result table exists. Real errors give other statuses (2 for a
# missing index, 137 for a kill)
KMA_OK_STATUS = (0, 95)

logger = None   # RGI's own logger (app.BWT.logger), set in main()


def fail(message, status):
    """Stop rgi with a failed child's exit status, 128 + signal for a killed one."""
    sys.stderr.write(f'rgi_bwt_kma_fix.py: {message}\n')
    sys.exit(128 - status if status < 0 else status or 1)


def run_filtered(cmd, output_sam_file, awk_program):
    """Run an aligner command line (SAM on stdout) through an awk filter.

    Returns the aligner's and the filter's exit statuses.
    """
    with open(output_sam_file, 'w') as out:
        aligner = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        awk = subprocess.Popen(['awk', '-F', '\t', awk_program], stdin=aligner.stdout, stdout=out)
        aligner.stdout.close()   # awk holds the pipe now
        awk_status, aligner_status = awk.wait(), aligner.wait()
    return aligner_status, awk_status


def run_kma(cmd, output_sam_file):
    """Run a KMA command line and keep the placed reads, grouped by name.

    Success needs an accepted KMA exit status, a successful filter, and
    KMA's result table (<output_sam_file>.temp.res, written when KMA closes
    its files; rgi bwt reads it later for the read coverage).
    """
    kma_status, awk_status = run_filtered(cmd, output_sam_file, KMA_PLACED_READS_AWK)
    # checked first: a failed filter closes the pipe, and KMA then dies of SIGPIPE
    if awk_status != 0:
        fail(f'filtering the KMA SAM failed (awk exit status {awk_status})', awk_status)
    if kma_status not in KMA_OK_STATUS or not os.path.exists(f'{output_sam_file}.temp.res'):
        fail(f'KMA failed (exit status {kma_status}): {cmd}', kma_status)


def run_aligner(name, cmd, output_sam_file):
    """Run a bowtie2/bwa command line and keep the records with a reference name."""
    status, awk_status = run_filtered(cmd, output_sam_file, ALIGNED_READS_AWK)
    if awk_status != 0:
        fail(f'filtering the {name} SAM failed (awk exit status {awk_status})', awk_status)
    if status != 0:
        fail(f'{name} failed (exit status {status}): {cmd}', status)


def unlink_staged_index(index_directory):
    """Remove the task's symlinks to the index files <index_directory>.*."""
    directory, prefix = os.path.split(index_directory)
    for path in glob.glob(os.path.join(glob.escape(directory), glob.escape(prefix) + '.*')):
        if os.path.islink(path):
            os.unlink(path)


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


def align_bowtie2(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_bowtie2 (RGI 6.0.8): reads with a reference, exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    logger.info("align reads -1 {} -2 {} to {}".format(self.read_one, self.read_two, reference_genome))
    run_aligner('bowtie2', "bowtie2 --quiet --very-sensitive-local --threads {threads} -x {index_directory} "
                "-1 {read_one} -2 {read_two}".format(
                    threads=self.threads, index_directory=index_directory, read_one=self.read_one,
                    read_two=self.read_two), output_sam_file)


def align_bowtie2_unpaired(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_bowtie2_unpaired (RGI 6.0.8): reads with a reference, exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    run_aligner('bowtie2', "bowtie2 --very-sensitive-local --threads {threads} -x {index_directory} "
                "-U {unpaired_reads}".format(
                    threads=self.threads, index_directory=index_directory,
                    unpaired_reads=self.read_one), output_sam_file)


def align_bwa_paired_end_mapping(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_bwa_paired_end_mapping (RGI 6.0.8): reads with a reference, exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    run_aligner('bwa', "bwa mem -t {threads} {index_directory} {read_one} {read_two}".format(
        threads=self.threads, index_directory=index_directory, read_one=self.read_one,
        read_two=self.read_two), output_sam_file)


def align_bwa_single_end_mapping(self, reference_genome, index_directory, output_sam_file):
    """Fixed copy of BWT.align_bwa_single_end_mapping (RGI 6.0.8): reads with a reference, exit check."""
    self.check_index(index_directory=index_directory, reference_genome=reference_genome)
    run_aligner('bwa', "bwa mem -M -t {threads} {index_directory} {read_one}".format(
        threads=self.threads, index_directory=index_directory, read_one=self.read_one), output_sam_file)


REPLACEMENTS = (
    ('align_kma', align_kma),
    ('align_kma_interleaved', align_kma_interleaved),
    ('align_bowtie2', align_bowtie2),
    ('align_bowtie2_unpaired', align_bowtie2_unpaired),
    ('align_bwa_paired_end_mapping', align_bwa_paired_end_mapping),
    ('align_bwa_single_end_mapping', align_bwa_single_end_mapping),
)


def patch(bwt_cls):
    """Replace the BWT methods after checking they are the 6.0.8 ones."""
    for name in RGI_608_METHOD_SHA256:
        try:
            installed = inspect.getsource(getattr(bwt_cls, name))
        except (OSError, TypeError, AttributeError):
            installed = ''   # no source to verify: treat as not the 6.0.8 method
        digest = hashlib.sha256(installed.encode()).hexdigest()
        if digest != RGI_608_METHOD_SHA256[name]:
            sys.exit(
                f'rgi_bwt_kma_fix.py: the installed RGI BWT.{name} is not the '
                f'6.0.8 method this patch targets (sha256 {digest}). Check '
                'whether this RGI version still writes every read to the '
                'SAM, then update or delete bin/rgi_bwt_kma_fix.py.')
    for name, replacement in REPLACEMENTS:
        setattr(bwt_cls, name, replacement)

    rgi_create_index = bwt_cls.create_index

    def create_index(self, index_directory, reference_genome):
        """BWT.create_index (RGI 6.0.8) after unlinking the staged index files."""
        unlink_staged_index(index_directory)
        rgi_create_index(self, index_directory=index_directory, reference_genome=reference_genome)

    bwt_cls.create_index = create_index


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
