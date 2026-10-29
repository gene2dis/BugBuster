#!/usr/bin/env python3
"""Contig -> refined bin -> the bin's GTDB-Tk classification, for one sample.

usage: bin_contig_labels.py <sample> <bins dir> [gtdbtk summary ...]

Writes <sample>_contig_bins.tsv with columns contig, bin,
bin_gtdb_classification. contig is the first word of a FASTA header, bin the
bin FASTA's name without .fa (GTDB-Tk's user_genome), and the classification
is GTDB-Tk's string, untouched (empty when GTDB-Tk has no row for the bin). A
contig already placed in an earlier bin is skipped with a warning (refined
bins should never share contigs). A bins directory without .fa files (only a
SKIPPED/FAILED marker) gives a header-only table.
"""
import glob
import os
import sys


def gtdb_classes(paths):
    """{user_genome: classification} from GTDB-Tk bac120/ar53 summaries."""
    classes = {}
    for path in paths:
        with open(path) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            if 'user_genome' not in header or 'classification' not in header:
                sys.exit(f"ERROR: {path} is not a GTDB-Tk summary (no user_genome/classification columns)")
            genome_col, class_col = header.index('user_genome'), header.index('classification')
            for line in fh:
                fields = line.rstrip('\n').split('\t')
                if len(fields) > max(genome_col, class_col) and fields[genome_col]:
                    classes[fields[genome_col]] = fields[class_col]
    return classes


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    sample, bins_dir, gtdb_files = sys.argv[1], sys.argv[2], sys.argv[3:]
    classes = gtdb_classes(gtdb_files)
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
                rows.append((contig, bin_id, classes.get(bin_id, '')))
    with open(f'{sample}_contig_bins.tsv', 'w') as out:
        out.write('contig\tbin\tbin_gtdb_classification\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')
    print(f"{sample}: {len(rows)} contigs in {len({r[1] for r in rows})} bins")


if __name__ == '__main__':
    main()
