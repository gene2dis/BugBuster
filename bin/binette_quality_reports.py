#!/usr/bin/env python3
"""
Binette's CheckM2 reports in CHECKM2_BATCH's per-binner format.

Binette runs CheckM2 on every input bin and every bin it selects
(input_bins_quality_reports/*.tsv, final_bins_quality_reports.tsv), so when
bins are refined the pipeline takes quality from it instead of a separate
CheckM2 run (performance review M7). This script writes the files
CHECKM2_BATCH would have written, with CheckM2 1.1.0's header:

    <sample>_<binner>_quality_report.tsv   one per staged binner directory
    <sample>_binette_quality_report.tsv    the refined bins

Binette reports completeness, contamination, model, coding density, N50,
size and contig count. Gene counts and mean protein length come from
Binette's own gene calls (temporary_files/assembly_proteins.faa.gz, one
Pyrodigal run on the whole assembly, not CheckM2's per-bin Prodigal);
GC content and the longest contig from the bin FASTA. Translation table and
notes are left empty. A binner without bins, or a sample Binette did not run
on, gets a header-only report.

The binner directories are passed under their binner's name (metabat,
semibin, comebin), as BINETTE stages them for Binette. Binette merges bins
with the same contigs across binners and lists the merged bin under one
binner's name in every binner's report, so each row is matched to the
binner's own bin by contig set.

Usage:
    binette_quality_reports.py --sample ID --final-bins DIR \
        [--binette-out DIR] BINNER_DIR [BINNER_DIR ...]
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

HEADER = [
    'Name', 'Completeness', 'Contamination', 'Completeness_Model_Used',
    'Translation_Table_Used', 'Coding_Density', 'Contig_N50',
    'Average_Gene_Length', 'Genome_Size', 'GC_Content',
    'Total_Coding_Sequences', 'Total_Contigs', 'Max_Contig_Length',
    'Additional_Notes',
]
FASTA_EXTENSIONS = ('.fa', '.fasta', '.fna')
# BINETTE stages each binner's directory under the binner's name (the name
# bin/identify_bin_dir.sh reports)
BINNERS = ('metabat', 'semibin', 'comebin')


def binner_of(dir_name):
    if dir_name not in BINNERS:
        sys.exit(f"ERROR: unrecognized binner directory '{dir_name}' (expected one of {', '.join(BINNERS)})")
    return dir_name


def fasta_bins(bin_dir):
    """Bin name (file name without extension) -> FASTA path."""
    return {
        f.name.rsplit('.', 1)[0]: f
        for f in sorted(Path(bin_dir).iterdir())
        if f.name.endswith(FASTA_EXTENSIONS)
    }


def fasta_stats(path):
    """Contig ids (first header word), GC fraction and longest contig."""
    contigs, gc, total, longest, length = [], 0, 0, 0, 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                contigs.append(line[1:].split()[0])
                longest, length = max(longest, length), 0
            else:
                seq = line.strip().upper()
                gc += seq.count('G') + seq.count('C')
                total += len(seq)
                length += len(seq)
    return contigs, (gc / total if total else 0.0), max(longest, length)


def gene_stats(proteins):
    """Contig id -> [genes, total protein length] from a Prodigal-format FAA
    (>contig_geneN # start # end ...)."""
    stats = {}
    if not proteins or not Path(proteins).exists():
        return None
    contig, length = None, 0

    def flush():
        if contig is not None:
            entry = stats.setdefault(contig, [0, 0])
            entry[0] += 1
            entry[1] += length

    with gzip.open(proteins, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                flush()
                contig, length = line[1:].split()[0].rsplit('_', 1)[0], 0
            else:
                length += len(line.strip().rstrip('*'))
    flush()
    return stats


def read_binette(report):
    with open(report) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def contig_set(fasta):
    with open(fasta) as f:
        return frozenset(line[1:].split()[0] for line in f if line.startswith('>'))


def checkm2_row(binette_row, fasta, genes, name=None):
    contigs, gc, longest = fasta_stats(fasta)
    n_genes = total_len = ''
    if genes is not None:
        counts = [genes.get(c, [0, 0]) for c in contigs]
        n = sum(c[0] for c in counts)
        n_genes = str(n)
        total_len = str(sum(c[1] for c in counts) / n) if n else '0'
    return [
        name or binette_row['name'], binette_row['completeness'], binette_row['contamination'],
        binette_row['checkm2_model'], '', binette_row['coding_density'], binette_row['N50'],
        total_len, binette_row['size'], f"{gc:.2f}", n_genes,
        binette_row['contig_count'], str(longest), '',
    ]


def write_report(path, rows):
    with open(path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(HEADER)
        writer.writerows(sorted(rows, key=lambda r: r[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--sample', required=True)
    parser.add_argument('--binette-out', help="Binette's output directory (absent when it did not run)")
    parser.add_argument('--final-bins', required=True, help='Directory of the refined bins')
    parser.add_argument('bin_dirs', nargs='+', help='Staged binner directories')
    args = parser.parse_args()

    out_dir = Path(args.binette_out) if args.binette_out else None
    ran = out_dir is not None and (out_dir / 'final_bins_quality_reports.tsv').exists()
    genes = gene_stats(out_dir / 'temporary_files' / 'assembly_proteins.faa.gz') if ran else None
    if ran and genes is None:
        print("WARNING: Binette's protein file is missing; gene counts left empty", file=sys.stderr)

    input_reports = {}
    if ran:
        # input_bins_<n>.<set name>.tsv; a single set is named after its path
        # (sets/metabat -> sets_metabat)
        for report in (out_dir / 'input_bins_quality_reports').glob('input_bins_*.tsv'):
            set_name = report.name[:-4].split('.', 1)[1]
            input_reports[set_name.split('_')[-1]] = report

    # Every binner's bins, to find the FASTA of a bin Binette reports under
    # another binner's name
    all_fastas = {}
    for bin_dir in args.bin_dirs:
        for bin_name, fasta in fasta_bins(bin_dir).items():
            all_fastas.setdefault(bin_name, fasta)

    for bin_dir in args.bin_dirs:
        binner = binner_of(Path(bin_dir).name)
        fastas = fasta_bins(bin_dir)
        rows = []
        if ran and fastas:
            if binner not in input_reports:
                sys.exit(f"ERROR: Binette wrote no report for the {binner} bins in {bin_dir}")
            own = {contig_set(fasta): name for name, fasta in fastas.items()}
            for row in read_binette(input_reports[binner]):
                fasta = fastas.get(row['name']) or all_fastas.get(row['name'])
                if fasta is None:
                    sys.exit(f"ERROR: Binette reports bin '{row['name']}' that is not in any binner directory")
                name = own.get(contig_set(fasta))
                if name is None:
                    sys.exit(f"ERROR: Binette reports bin '{row['name']}' for {binner}, but no {binner} bin has its contigs")
                rows.append(checkm2_row(row, fastas[name], genes, name))
        write_report(f"{args.sample}_{binner}_quality_report.tsv", rows)

    rows = []
    if ran:
        fastas = fasta_bins(args.final_bins)
        for row in read_binette(out_dir / 'final_bins_quality_reports.tsv'):
            if row['name'] not in fastas:
                sys.exit(f"ERROR: Binette reports final bin '{row['name']}' that is not in {args.final_bins}")
            rows.append(checkm2_row(row, fastas[row['name']], genes))
    write_report(f"{args.sample}_binette_quality_report.tsv", rows)


if __name__ == '__main__':
    main()
