#!/usr/bin/env python3
"""
Bin Taxonomy Report Generator
Creates bar plots and summary tables for MAG taxonomy at phylum level.
Python migration of Bins_tax.R
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style as fs  # noqa: E402

fs.setup()
import matplotlib.pyplot as plt  # noqa: E402


# Columns of MAGs_tax_summary.csv
SUMMARY_COLUMNS = [
    'sample', 'bin_name', 'closest_genome_reference',
    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
]


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
    """Main function to generate bin taxonomy report."""
    
    # Find GTDB-Tk taxonomy files
    tax_files = glob.glob('*gtdbtk*.tsv')
    if not tax_files:
        print("ERROR: No GTDB-Tk taxonomy files found (*gtdbtk*.tsv)", file=sys.stderr)
        sys.exit(1)
    
    # Read and process all taxonomy files
    all_tax_reports = []
    for file in tax_files:
        sample_name = file.split('_gtdbtk_')[0]
        df = pd.read_csv(file, sep='\t')
        
        # Skip empty files (only header, no data rows)
        if len(df) == 0:
            print(f"WARNING: Skipping empty file {file} (no bins classified)")
            continue
        
        # Parse taxonomy classification
        taxonomy_parsed = df['classification'].apply(parse_gtdb_taxonomy)
        taxonomy_df = pd.DataFrame(taxonomy_parsed.tolist())
        
        # Use correct column name from GTDB-Tk output
        # GTDB-Tk uses 'fastani_reference' not 'closest_genome_reference'
        reference_col = 'fastani_reference' if 'fastani_reference' in df.columns else 'closest_genome_reference'
        
        # Combine with original data
        result = pd.concat([
            df[['user_genome', 'classification', reference_col]], 
            taxonomy_df
        ], axis=1)
        result.rename(columns={reference_col: 'closest_genome_reference'}, inplace=True)
        result['sample'] = sample_name
        
        all_tax_reports.append(result)
    
    # No sample has a classified bin: binning found nothing usable on this
    # data (SemiBin and Binette now finish with empty bin sets). Warn and
    # write a header-only summary and a placeholder plot so the run can
    # finish; a missing GTDB-Tk file above still fails
    if not all_tax_reports:
        print(
            "WARNING: all GTDB-Tk summaries are empty (header-only) - no bins "
            "were classified. Upstream binning produced no refined bins for "
            "any sample; writing an empty taxonomy summary.",
            file=sys.stderr
        )
        pd.DataFrame(columns=SUMMARY_COLUMNS).to_csv('MAGs_tax_summary.csv', index=False)
        fs.placeholder('.', 'MAGs_tax_plot', 'No MAGs were classified')
        print("Generated an empty MAGs_tax_summary.csv and a placeholder MAGs_tax_plot")
        return

    # Combine all reports
    all_tax_bin_reports_full = pd.concat(all_tax_reports, ignore_index=True)
    
    # Select final columns and rename
    all_tax_bin_reports_csv = all_tax_bin_reports_full[[
        'sample', 'user_genome', 'closest_genome_reference',
        'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
    ]].rename(columns={'user_genome': 'bin_name'})
    assert list(all_tax_bin_reports_csv.columns) == SUMMARY_COLUMNS
    
    # Write CSV output
    all_tax_bin_reports_csv.to_csv('MAGs_tax_summary.csv', index=False)

    # Bins per GTDB phylum; a bin without a phylum is counted as unclassified, not dropped
    phyla = all_tax_bin_reports_csv['Phylum'].fillna('unclassified').replace('', 'unclassified')
    phylum_counts = phyla.value_counts().rename_axis('Phylum').reset_index(name='Bin count')
    phylum_counts = phylum_counts.sort_values(['Bin count', 'Phylum'], ascending=[False, True], kind='stable')
    phylum_counts = pd.concat([phylum_counts[phylum_counts['Phylum'] != 'unclassified'],
                               phylum_counts[phylum_counts['Phylum'] == 'unclassified']], ignore_index=True)
    n_phyla = len(phylum_counts)
    fig, ax = plt.subplots(figsize=fs.figsize('single', n_phyla * 3.2 + 16), layout='constrained')
    y = np.arange(n_phyla)
    colors = [fs.OTHER if p == 'unclassified' else fs.CATEGORY[0] for p in phylum_counts['Phylum']]
    ax.barh(y, phylum_counts['Bin count'], color=colors, height=0.7, lw=0)
    for i, v in enumerate(phylum_counts['Bin count']):
        ax.text(v, i, f" {v}", va='center', fontsize=fs.SMALL_PT)
    ax.set_yticks(y, phylum_counts['Phylum'])
    ax.set_ylim(n_phyla - 0.5, -0.5)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('MAGs (GTDB-Tk phylum)')
    ax.set_xlim(0, phylum_counts['Bin count'].max() * 1.15)
    fs.save(fig, '.', 'MAGs_tax_plot', phylum_counts.rename(columns={'Bin count': 'mags'}))
    print("Successfully generated MAGs_tax_plot")

    print(f"Successfully generated MAGs_tax_summary.csv with {len(all_tax_bin_reports_csv)} MAGs")
    print(f"  - Unique phyla: {n_phyla}")
    print(f"  - Phylum distribution: {dict(zip(phylum_counts['Phylum'], phylum_counts['Bin count']))}")


if __name__ == '__main__':
    main()
