#!/usr/bin/env python3
"""Contig -> refined bin, for one sample.

usage: bin_contig_labels.py <sample> <bins dir>

Writes <sample>_contig_bins.tsv with columns contig, bin. contig is the first
word of a FASTA header and bin the bin FASTA's name without .fa (GTDB-Tk's
user_genome, which ARG_CONTIG_LEVEL_REPORT joins for the bin's
classification). A contig already placed in an earlier bin is skipped with a
warning (refined bins should never share contigs). A bins directory without
.fa files (only a SKIPPED/FAILED marker) gives a header-only table. Needs only
the bins, so the bin-level ARG split does not wait for GTDB-Tk.
"""
import glob
import os
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sample, bins_dir = sys.argv[1], sys.argv[2]
    seen = {}
    rows = []
    for fasta in sorted(glob.glob(os.path.join(bins_dir, '*.fa'))):
        bin_id = os.path.basename(fasta)[:-len('.fa')]
        with open(fasta) as fh:
            for line in fh:
                if not line.startswith('>'):
                    continue
                contig = line[1:].split()[0]
                if contig in seen:
                    print(f"WARNING: contig {contig} of bin {bin_id} is already in bin {seen[contig]} - skipped",
                          file=sys.stderr)
                    continue
                seen[contig] = bin_id
                rows.append((contig, bin_id))
    with open(f'{sample}_contig_bins.tsv', 'w') as out:
        out.write('contig\tbin\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')
    print(f"{sample}: {len(rows)} contigs in {len({r[1] for r in rows})} bins")


if __name__ == '__main__':
    main()
