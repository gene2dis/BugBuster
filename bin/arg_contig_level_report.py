#!/usr/bin/env python3
"""
ARG Contig Level Report Generator
Joins, per contig, the BlobTools table (length, GC, coverage), Metabuli's
classification (seven ranks from its lineage, the rank it called at and its
score) and the DeepARG predictions. With binning, the contig's refined bin
and that bin's GTDB-Tk classification (<sample>_contig_bins.tsv) come as two
more columns, next to the NCBI ranks and never filling them. With geNomad
(--genomad_prediction), each contig gets its mobile class (plasmid or virus
call of geNomad's summaries, else chromosome), its aggregated plasmid and
virus scores, virus taxonomy and plasmid conjugation genes, and each ARG gene
the provirus regions it overlaps (gene coordinates from the Pyrodigal GFF).
Python migration of Contig_arg_unify.R
"""

import pandas as pd
import glob
import gzip
import re
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


GENOMAD_COLS = ['mobile_class', 'plasmid_score', 'virus_score', 'provirus', 'virus_taxonomy',
                'plasmid_conjugation_genes']
# geNomad names its outputs after the contigs file, <sample>_filtered_contigs
GENOMAD_SUFFIX = '_filtered_contigs_'


def read_tsv(path):
    """A TSV as strings; an empty file (a stub, or a table without rows) gives an empty frame."""
    try:
        return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def genomad_tables():
    """{sample: (contig-level frame, {contig: [(start, end, taxonomy)] of its proviruses})}."""
    out = {}
    for path in glob.glob(f'*{GENOMAD_SUFFIX}plasmid_summary.tsv'):
        sample = path[:-len(f'{GENOMAD_SUFFIX}plasmid_summary.tsv')]
        plasmids = read_tsv(path)
        viruses = read_tsv(f'{sample}{GENOMAD_SUFFIX}virus_summary.tsv')
        scores = read_tsv(f'{sample}{GENOMAD_SUFFIX}aggregated_classification.tsv')
        contig = pd.DataFrame(columns=['contig_name'])
        if not scores.empty:
            contig = scores[['seq_name', 'plasmid_score', 'virus_score']].rename(columns={'seq_name': 'contig_name'})
        proviruses = {}
        whole_viruses = pd.DataFrame(columns=['contig_name', 'virus_taxonomy'])
        if not viruses.empty:
            is_provirus = viruses['seq_name'].str.contains('|provirus_', regex=False)
            for name, taxonomy in viruses.loc[is_provirus, ['seq_name', 'taxonomy']].itertuples(index=False):
                host, region = name.split('|', 1)
                _, start, end = region.split('_')
                proviruses.setdefault(host, []).append((int(start), int(end), taxonomy))
            whole_viruses = viruses.loc[~is_provirus, ['seq_name', 'taxonomy']].rename(
                columns={'seq_name': 'contig_name', 'taxonomy': 'virus_taxonomy'})
        plasmid_names = set(plasmids['seq_name']) if not plasmids.empty else set()
        contig = contig.merge(whole_viruses, on='contig_name', how='outer')
        if not plasmids.empty:
            contig = contig.merge(plasmids[['seq_name', 'conjugation_genes']].rename(
                columns={'seq_name': 'contig_name', 'conjugation_genes': 'plasmid_conjugation_genes'}),
                on='contig_name', how='outer')
        out[sample] = (contig, proviruses, plasmid_names, set(whole_viruses['contig_name']))
    return out


def gene_coordinates(sample):
    """{gene id: (start, end)} from the sample's Pyrodigal GFF (<sample>.gff.gz)."""
    coords = {}
    path = f'{sample}.gff.gz'
    try:
        with gzip.open(path, 'rt') as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if line.startswith('#') or len(f) < 9 or f[2] != 'CDS':
                    continue
                gene = re.search(r'ID=([^;]+)', f[8])
                if gene:
                    coords[gene.group(1)] = (int(f[3]), int(f[4]))
    except FileNotFoundError:
        sys.exit(f"ERROR: geNomad tables for {sample} but no Pyrodigal GFF ({path}) for the ARG gene coordinates")
    return coords


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

    # Bin membership and the bin's GTDB-Tk classification (binning on)
    bin_files = glob.glob('*_contig_bins.tsv')
    bin_cols = []
    if bin_files:
        bins = pd.concat(
            [pd.read_csv(f, sep='\t', dtype=str, keep_default_na=False).assign(Sample=f[:-len('_contig_bins.tsv')])
             for f in bin_files],
            ignore_index=True).rename(columns={'contig': 'contig_name'})
        all_blob_reports = all_blob_reports.merge(
            bins[['Sample', 'contig_name', 'bin', 'bin_gtdb_classification']],
            on=['Sample', 'contig_name'], how='left')
        bin_cols = ['bin', 'bin_gtdb_classification']

    # geNomad: contig-level calls and scores (--genomad_prediction)
    genomad = genomad_tables()
    genomad_cols = []
    if genomad:
        per_sample = []
        for sample, (contig, _proviruses, plasmid_names, virus_names) in genomad.items():
            part = all_blob_reports.loc[all_blob_reports['Sample'] == sample, ['Sample', 'contig_name']].merge(
                contig, on='contig_name', how='left')
            part['mobile_class'] = 'chromosome'
            part.loc[part['contig_name'].isin(plasmid_names), 'mobile_class'] = 'plasmid'
            part.loc[part['contig_name'].isin(virus_names), 'mobile_class'] = 'virus'
            per_sample.append(part)
        all_blob_reports = all_blob_reports.merge(
            pd.concat(per_sample, ignore_index=True).reindex(
                columns=['Sample', 'contig_name', 'mobile_class', 'plasmid_score', 'virus_score',
                         'virus_taxonomy', 'plasmid_conjugation_genes']),
            on=['Sample', 'contig_name'], how='left')
        genomad_cols = GENOMAD_COLS

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
        'Sample', 'contig_name', 'read_id', '#ARG', 'predicted_ARG-class',
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

    # geNomad proviruses overlapping each ARG gene: their coordinates, and
    # their taxonomy when the contig itself is not a virus call
    if genomad:
        all_data['provirus'] = ''
        for sample, (_contig, proviruses, _plasmids, _viruses) in genomad.items():
            coords = gene_coordinates(sample)
            rows = all_data.index[(all_data['Sample'] == sample) & all_data['read_id'].notna()
                                  & all_data['contig_name'].isin(list(proviruses))]
            for i in rows:
                gene = coords.get(all_data.at[i, 'read_id'])
                if gene is None:
                    continue
                hits = [(a, b, tax) for a, b, tax in proviruses[all_data.at[i, 'contig_name']]
                        if gene[0] <= b and gene[1] >= a]
                if hits:
                    all_data.at[i, 'provirus'] = ';'.join(f'{a}-{b}' for a, b, _ in hits)
                    current = all_data.at[i, 'virus_taxonomy']
                    if pd.isna(current) or current == '':
                        all_data.at[i, 'virus_taxonomy'] = ' | '.join(tax for _, _, tax in hits)

    # Reorder columns
    all_data = all_data[[
        'Sample', 'contig_name', 'contig_length', 'GC', 'coverage'
    ] + RANKS + [
        'taxonomy_rank', 'taxonomy_score'
    ] + bin_cols + genomad_cols + [
        'ARG_gene', 'predicted_ARG_class', 'best_hit', 'probability',
        'identity', 'alignment_evalue', 'counts'
    ]]

    # Write output (tab-separated, as the file name says)
    all_data.to_csv('Contig_tax_and_arg_prediction.tsv', sep='\t', index=False)

    print(f"Successfully generated Contig_tax_and_arg_prediction.tsv with {len(all_data)} contigs")
    print(f"  - Contigs with ARG predictions: {all_data['ARG_gene'].notna().sum()}")


if __name__ == '__main__':
    main()
