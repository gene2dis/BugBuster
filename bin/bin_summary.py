#!/usr/bin/env python3
"""
Bin Summary Report Generator
Combines bin quality, taxonomy, and depth information into a single summary table.
Python migration of Bin_summary.R
"""

import pandas as pd
import glob
import sys
from pathlib import Path


def parse_gtdb_taxonomy(classification):
    """Parse GTDB taxonomy string into individual ranks."""
    if pd.isna(classification):
        return {
            'Domain': None, 'Phylum': None, 'Class': None,
            'Order': None, 'Family': None, 'Genus': None, 'Species': None
        }
    
    ranks = classification.split(';')
    taxonomy = {}
    rank_names = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    for i, rank_name in enumerate(rank_names):
        if i < len(ranks):
            # Extract value after '__' prefix (e.g., 'd__Bacteria' -> 'Bacteria')
            parts = ranks[i].split('__')
            taxonomy[rank_name] = parts[1] if len(parts) > 1 and parts[1] else None
        else:
            taxonomy[rank_name] = None
    
    return taxonomy


def main():
    """Main function to generate bin summary report."""
    
    # Read bin depth files
    depth_files = glob.glob('*bin_depth.tsv')
    if not depth_files:
        print("ERROR: No bin depth files found (*bin_depth.tsv)", file=sys.stderr)
        sys.exit(1)
    
    depth_dfs = []
    for file in depth_files:
        df = pd.read_csv(file, sep='\t')
        depth_dfs.append(df)
    
    bin_depth = pd.concat(depth_dfs, ignore_index=True)
    
    # Add suffix to sample names and pivot
    bin_depth['sample'] = bin_depth['sample'] + '_avgcov'
    bin_depth_wide = bin_depth.pivot_table(
        index=['bin_id', 'Total_contigs'],
        columns='sample',
        values='average_cov',
        aggfunc='first'
    ).reset_index()
    
    # Separate Total_contigs for later join
    bin_tmp = bin_depth_wide[['bin_id', 'Total_contigs']].copy()
    
    # Read ALL bin quality reports (per-binner, per-sample); previously only the
    # first glob hit was read, dropping every other staged report (audit #11)
    quality_files = glob.glob('*quality_report.tsv')
    if not quality_files:
        print("ERROR: No quality report files found (*quality_report.tsv)", file=sys.stderr)
        sys.exit(1)

    quality_dfs = []
    for file in sorted(quality_files):
        try:
            quality_dfs.append(pd.read_csv(file, sep='\t'))
        except pd.errors.EmptyDataError:
            print(f"WARNING: skipping empty quality report {file}", file=sys.stderr)

    bin_quality = pd.concat(quality_dfs, ignore_index=True) if quality_dfs else pd.DataFrame()
    if bin_quality.empty:
        print("ERROR: all quality reports are empty - no bins were assessed", file=sys.stderr)
        sys.exit(1)

    bin_quality = bin_quality[[
        'Name', 'Completeness', 'Contamination',
        'Contig_N50', 'Genome_Size', 'Total_Coding_Sequences'
    ]]

    # Read ALL bin taxonomy summaries (bac120 + ar53); GTDB-Tk writes a zero-byte
    # or header-only ar53 file when no archaea are found
    tax_files = glob.glob('*gtdbtk*.tsv')
    if not tax_files:
        print("ERROR: No GTDB-Tk taxonomy files found (*gtdbtk*.tsv)", file=sys.stderr)
        sys.exit(1)

    tax_dfs = []
    for file in sorted(tax_files):
        try:
            tax_dfs.append(pd.read_csv(file, sep='\t'))
        except pd.errors.EmptyDataError:
            print(f"WARNING: skipping empty taxonomy file {file}", file=sys.stderr)

    bin_tax = pd.concat(tax_dfs, ignore_index=True) if tax_dfs else pd.DataFrame()
    if bin_tax.empty:
        print("ERROR: all GTDB-Tk summaries are empty - no bins were classified", file=sys.stderr)
        sys.exit(1)
    
    # Parse taxonomy classification
    taxonomy_parsed = bin_tax['classification'].apply(parse_gtdb_taxonomy)
    taxonomy_df = pd.DataFrame(taxonomy_parsed.tolist())
    
    bin_tax = pd.concat([bin_tax[['user_genome']], taxonomy_df], axis=1)
    
    # Fill missing species with "Genus sp."
    bin_tax['Species'] = bin_tax.apply(
        lambda row: f"{row['Genus']} sp." if pd.isna(row['Species']) and pd.notna(row['Genus']) else row['Species'],
        axis=1
    )
    
    # Join all dataframes
    bin_reports = bin_quality.merge(
        bin_tmp,
        left_on='Name',
        right_on='bin_id',
        how='left'
    ).drop(columns=['bin_id'])
    
    bin_reports = bin_reports.merge(
        bin_tax,
        left_on='Name',
        right_on='user_genome',
        how='left'
    ).drop(columns=['user_genome'])
    
    # Join with depth data (excluding Total_contigs which is already present)
    depth_cols = [col for col in bin_depth_wide.columns if col not in ['bin_id', 'Total_contigs']]
    bin_depth_for_join = bin_depth_wide[['bin_id'] + depth_cols]
    
    bin_reports = bin_reports.merge(
        bin_depth_for_join,
        left_on='Name',
        right_on='bin_id',
        how='left'
    ).drop(columns=['bin_id'])
    
    # Write output
    bin_reports.to_csv('Bin_summary.csv', index=False)
    print(f"Successfully generated Bin_summary.csv with {len(bin_reports)} bins")


if __name__ == '__main__':
    main()
