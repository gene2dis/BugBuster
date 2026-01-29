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
    
    # Read bin quality report
    quality_files = glob.glob('*quality_report.tsv')
    if not quality_files:
        print("ERROR: No quality report files found (*quality_report.tsv)", file=sys.stderr)
        sys.exit(1)
    
    bin_quality = pd.read_csv(quality_files[0], sep='\t')
    bin_quality = bin_quality[[
        'Name', 'Completeness', 'Contamination', 
        'Contig_N50', 'Genome_Size', 'Total_Coding_Sequences'
    ]]
    
    # Read bin taxonomy
    tax_files = glob.glob('*gtdbtk*.tsv')
    if not tax_files:
        print("ERROR: No GTDB-Tk taxonomy files found (*gtdbtk*.tsv)", file=sys.stderr)
        sys.exit(1)
    
    bin_tax = pd.read_csv(tax_files[0], sep='\t')
    
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
