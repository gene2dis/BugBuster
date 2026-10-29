#!/usr/bin/env python3
"""
ARG Contig Level Report Generator
Joins, per contig, the BlobTools table (length, GC, coverage), Metabuli's
classification (seven ranks from its lineage, the rank it called at and its
score) and the DeepARG predictions. Python migration of Contig_arg_unify.R
"""

import pandas as pd
import glob
import sys

RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
# Metabuli lineage prefixes (--lineage 1); 'd' is NCBI's 'domain' (the
# superkingdom before 2025) or GTDB's domain
LINEAGE_PREFIX = {'d': 'superkingdom', 'p': 'phylum', 'c': 'class', 'o': 'order',
                  'f': 'family', 'g': 'genus', 's': 'species'}


def read_metabuli(path, sample):
    """Metabuli classify --lineage 1 output -> one row per contig.

    Parsed by position: lines end with a tab (the taxID:match_count field is
    empty above species), which a header-driven reader would turn into an
    extra column.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            row = {'Sample': sample, 'contig_name': f[1]}
            if f[0] == '1':
                for node in f[7].split(';'):
                    prefix, _, name = node.partition('_')
                    if prefix in LINEAGE_PREFIX:
                        row[LINEAGE_PREFIX[prefix]] = name
                    # NCBI's Viruses (10239) is an 'acellular root' since 2025,
                    # with no 'd_' node: it is the viruses' superkingdom
                    if name == 'Viruses':
                        row.setdefault('superkingdom', 'Viruses')
                row['taxonomy_rank'] = f[6]
                row['taxonomy_score'] = float(f[4])
            else:
                row['taxonomy_rank'] = 'unclassified'
                row['taxonomy_score'] = 0.0
            rows.append(row)
    return pd.DataFrame(rows, columns=['Sample', 'contig_name'] + RANKS + ['taxonomy_rank', 'taxonomy_score'])


def main():
    """Main function to generate contig-level ARG and taxonomy report."""

    # Read blob table files
    blob_files = glob.glob('*Blob_table*')
    if not blob_files:
        print("ERROR: No blob table files found (*Blob_table*)", file=sys.stderr)
        sys.exit(1)

    blob_dfs = []
    for file in blob_files:
        sample_name = file.split('_Blob_table')[0]
        # Skip first 10 lines (blobtools header)
        df = pd.read_csv(file, sep='\t', skiprows=10)
        df['Sample'] = sample_name
        blob_dfs.append(df)

    all_blob_reports = pd.concat(blob_dfs, ignore_index=True)

    # Contig length, GC and coverage from BlobTools; the ranks come from
    # Metabuli (BlobTools has no class rank)
    all_blob_reports = all_blob_reports[[
        'Sample', '# name', 'length', 'GC', 'bam0'
    ]].rename(columns={
        '# name': 'contig_name',
        'length': 'contig_length',
        'bam0': 'coverage'
    })

    metabuli_files = glob.glob('*_classifications.tsv')
    if not metabuli_files:
        print("ERROR: No Metabuli files found (*_classifications.tsv)", file=sys.stderr)
        sys.exit(1)
    all_tax = pd.concat(
        [read_metabuli(f, f[:-len('_classifications.tsv')]) for f in metabuli_files],
        ignore_index=True)
    all_blob_reports = all_blob_reports.merge(all_tax, on=['Sample', 'contig_name'], how='left')

    # Read DeepARG files
    deeparg_files = glob.glob('*_contigs_deep_arg.out.mapping.ARG')
    if not deeparg_files:
        print("ERROR: No DeepARG files found (*_contigs_deep_arg.out.mapping.ARG)", file=sys.stderr)
        sys.exit(1)

    arg_dfs = []
    for file in deeparg_files:
        sample_name = file.split('_contigs_deep_arg')[0]
        df = pd.read_csv(file, sep='\t')
        df['Sample'] = sample_name
        arg_dfs.append(df)

    all_arg_reports = pd.concat(arg_dfs, ignore_index=True)

    # read_id is the Pyrodigal gene id, <contig>_<gene number>
    all_arg_reports['contig_name'] = all_arg_reports['read_id'].str.rsplit('_', n=1).str[0]

    # Select and rename columns
    all_arg_reports = all_arg_reports[[
        'Sample', 'contig_name', '#ARG', 'predicted_ARG-class',
        'best-hit', 'probability', 'identity', 'alignment-evalue', 'counts'
    ]].rename(columns={
        '#ARG': 'ARG_gene',
        'predicted_ARG-class': 'predicted_ARG_class',
        'best-hit': 'best_hit',
        'alignment-evalue': 'alignment_evalue'
    })

    # Left join blob reports with ARG reports
    all_data = all_blob_reports.merge(
        all_arg_reports,
        on=['Sample', 'contig_name'],
        how='left'
    )

    # Reorder columns
    all_data = all_data[[
        'Sample', 'contig_name', 'contig_length', 'GC', 'coverage'
    ] + RANKS + [
        'taxonomy_rank', 'taxonomy_score',
        'ARG_gene', 'predicted_ARG_class', 'best_hit', 'probability',
        'identity', 'alignment_evalue', 'counts'
    ]]

    # Write output (tab-separated, as the file name says)
    all_data.to_csv('Contig_tax_and_arg_prediction.tsv', sep='\t', index=False)

    print(f"Successfully generated Contig_tax_and_arg_prediction.tsv with {len(all_data)} contigs")
    print(f"  - Contigs with ARG predictions: {all_data['ARG_gene'].notna().sum()}")


if __name__ == '__main__':
    main()
