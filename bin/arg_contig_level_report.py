#!/usr/bin/env python3
"""
ARG Contig Level Report Generator
Combines contig taxonomy (from blobtools) with ARG predictions (from DeepARG).
Python migration of Contig_arg_unify.R
"""

import pandas as pd
import glob
import sys


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
    
    # Select and rename columns
    all_blob_reports = all_blob_reports[[
        'Sample', '# name', 'length', 'GC', 'bam0',
        'superkingdom.t.6%s', 'phylum.t.9%s', 'order.t.12%s',
        'family.t.15%s', 'genus.t.18%s', 'species.t.21%s'
    ]].rename(columns={
        '# name': 'contig_name',
        'length': 'contig_length',
        'bam0': 'coverage',
        'superkingdom.t.6%s': 'superkingdom',
        'phylum.t.9%s': 'phylum',
        'order.t.12%s': 'order',
        'family.t.15%s': 'family',
        'genus.t.18%s': 'genus',
        'species.t.21%s': 'species'
    })
    
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
    
    # Parse read_id to extract contig_name
    # Format: contig_1_2 -> extract first two parts
    all_arg_reports['contig_parts'] = all_arg_reports['read_id'].str.split('_')
    all_arg_reports['new_read_id_1'] = all_arg_reports['contig_parts'].str[0]
    all_arg_reports['new_read_id_2'] = all_arg_reports['contig_parts'].str[1]
    all_arg_reports['contig_name'] = (
        all_arg_reports['new_read_id_1'] + '_' + all_arg_reports['new_read_id_2']
    )
    
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
        'Sample', 'contig_name', 'contig_length', 'GC', 'coverage',
        'superkingdom', 'phylum', 'order', 'family', 'genus', 'species',
        'ARG_gene', 'predicted_ARG_class', 'best_hit', 'probability',
        'identity', 'alignment_evalue', 'counts'
    ]]
    
    # Write output
    all_data.to_csv('Contig_tax_and_arg_prediction.tsv', index=False)
    
    print(f"Successfully generated Contig_tax_and_arg_prediction.tsv with {len(all_data)} contigs")
    print(f"  - Contigs with ARG predictions: {all_data['ARG_gene'].notna().sum()}")


if __name__ == '__main__':
    main()
