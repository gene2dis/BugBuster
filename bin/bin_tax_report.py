#!/usr/bin/env python3
"""
Bin Taxonomy Report Generator
Creates bar plots and summary tables for MAG taxonomy at phylum level.
Python migration of Bins_tax.R
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import glob
import sys
import numpy as np


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


def generate_colors(n):
    """Generate a pastel color palette."""
    # Use a colormap similar to rcartocolor Pastel
    cmap = plt.cm.get_cmap('Pastel1')
    if n <= 9:
        colors = [cmap(i) for i in np.linspace(0, 1, 9)][:n]
    else:
        # For more colors, use Set3 which has more variety
        cmap = plt.cm.get_cmap('Set3')
        colors = [cmap(i) for i in np.linspace(0, 1, min(n, 12))]
        if n > 12:
            # Generate additional colors
            cmap2 = plt.cm.get_cmap('Pastel2')
            colors.extend([cmap2(i) for i in np.linspace(0, 1, n - 12)])
    return colors[:n]


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
    
    # Check if we have any data
    if not all_tax_reports:
        print("WARNING: No bins were classified by GTDB-Tk. Creating empty output files.")
        # Create empty output files
        empty_df = pd.DataFrame(columns=[
            'sample', 'bin_name', 'closest_genome_reference',
            'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
        ])
        empty_df.to_csv('MAGs_tax_summary.csv', index=False)
        
        # Create empty plot
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.text(0.5, 0.5, 'No MAGs classified', 
                ha='center', va='center', fontsize=20, transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        plt.savefig('MAGs_tax_plot.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print("Created empty output files (no bins to classify)")
        return
    
    # Combine all reports
    all_tax_bin_reports_full = pd.concat(all_tax_reports, ignore_index=True)
    
    # Select final columns and rename
    all_tax_bin_reports_csv = all_tax_bin_reports_full[[
        'sample', 'user_genome', 'closest_genome_reference',
        'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
    ]].rename(columns={'user_genome': 'bin_name'})
    
    # Count bins by phylum
    phylum_counts = all_tax_bin_reports_csv.groupby('Phylum').size().reset_index(name='Bin count')
    phylum_counts = phylum_counts.sort_values('Bin count', ascending=False)
    
    # Write CSV output
    all_tax_bin_reports_csv.to_csv('MAGs_tax_summary.csv', index=False)
    
    # Create the plot
    try:
        fig, ax = plt.subplots(figsize=(16, 8))
        
        # Generate colors
        n_phyla = len(phylum_counts)
        colors = generate_colors(n_phyla)
        color_map = dict(zip(phylum_counts['Phylum'], colors))
        
        # Create bar plot
        bars = ax.bar(range(len(phylum_counts)), phylum_counts['Bin count'],
                      color=[color_map[p] for p in phylum_counts['Phylum']],
                      edgecolor='black', linewidth=1)
        
        # Add value labels on top of bars
        for i, (idx, row) in enumerate(phylum_counts.iterrows()):
            ax.text(i, row['Bin count'], str(row['Bin count']),
                   ha='center', va='bottom', fontsize=20)
        
        # Set labels and formatting
        ax.set_xticks(range(len(phylum_counts)))
        ax.set_xticklabels(phylum_counts['Phylum'], rotation=40, ha='right', fontsize=20)
        ax.set_ylabel('MAGs Count', fontsize=25)
        ax.set_xlabel(None)
        
        # Style the plot
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.7, color='gray')
        ax.set_axisbelow(True)
        ax.tick_params(axis='y', labelsize=18)
        
        # Set background
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        
        plt.tight_layout()
        plt.savefig('MAGs_tax_plot.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print("Successfully generated MAGs_tax_plot.png")
        
    except Exception as e:
        print(f"ERROR generating plot: {e}", file=sys.stderr)
    
    print(f"Successfully generated MAGs_tax_summary.csv with {len(all_tax_bin_reports_csv)} MAGs")
    print(f"  - Unique phyla: {n_phyla}")
    print(f"  - Phylum distribution: {dict(zip(phylum_counts['Phylum'], phylum_counts['Bin count']))}")


if __name__ == '__main__':
    main()
