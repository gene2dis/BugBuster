#!/usr/bin/env python3
"""
Phylum Blob Plot Generator
Creates multi-panel blob plots showing contig GC vs coverage colored by phylum.
Python migration of Blobplot.R
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import numpy as np
import glob
import sys
import pickle


def generate_pastel_colors(n):
    """Generate pastel color palette similar to rcartocolor Pastel."""
    # Pastel colors
    pastel_colors = [
        '#66C5CC', '#F6CF71', '#F89C74', '#DCB0F2', '#87C55F',
        '#9EB9F3', '#FE88B1', '#C9DB74', '#8BE0A4', '#B497E7',
        '#D3B484', '#B3B3B3'
    ]
    if n <= len(pastel_colors):
        return pastel_colors[:n]
    else:
        # If more colors needed, cycle through
        return [pastel_colors[i % len(pastel_colors)] for i in range(n)]


def main():
    """Main function to generate phylum blob plot."""
    
    # Find blob table files
    blob_files = glob.glob('*Blob_table*')
    if not blob_files:
        print("ERROR: No blob table files found (*Blob_table*)", file=sys.stderr)
        sys.exit(1)
    
    # Read and process all blob files
    all_blob_reports = []
    for file in blob_files:
        sample_name = file.split('_Blob_table')[0]
        # Skip first 10 lines (blobtools header)
        df = pd.read_csv(file, sep='\t', skiprows=10)
        df['Sample'] = sample_name
        all_blob_reports.append(df)
    
    all_blob_reports = pd.concat(all_blob_reports, ignore_index=True)
    
    # Select and rename columns
    data = all_blob_reports[[
        '# name', 'length', 'GC', 'bam0', 'phylum.t.9%s'
    ]].rename(columns={
        '# name': 'contig_name',
        'bam0': 'coverage',
        'phylum.t.9%s': 'phylum'
    })
    
    # Calculate total length by phylum
    phylum_summary = data.groupby('phylum')['length'].sum().reset_index()
    phylum_summary.columns = ['phylum', 'length_sum']
    phylum_summary = phylum_summary.sort_values('length_sum', ascending=False)
    
    # Determine number of top phyla to show (max 9)
    num_top = min(9, len(phylum_summary))
    
    # Join and classify phyla
    data = data.merge(phylum_summary, on='phylum', how='left')
    threshold = phylum_summary.iloc[num_top - 1]['length_sum'] if num_top > 0 else 0
    data['phylum'] = data.apply(
        lambda row: 'other' if row['length_sum'] <= threshold else row['phylum'],
        axis=1
    )
    
    # Generate colors
    unique_phyla = sorted(data['phylum'].unique())
    colors_list = generate_pastel_colors(len(unique_phyla))
    colors = dict(zip(unique_phyla, colors_list))
    
    # Special colors for specific phyla
    if 'other' in colors:
        colors['other'] = 'gray'
    if 'no-hit' in colors:
        colors['no-hit'] = 'white'
    
    # Calculate percentages
    phylum_counts = data.groupby('phylum').size().reset_index(name='n')
    phylum_counts['percent'] = (phylum_counts['n'] / phylum_counts['n'].sum() * 100).round(2)
    
    # Create figure with custom layout
    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.15)
    
    # Top panel: Blob plot
    ax_blob = fig.add_subplot(gs[0])
    
    # Log transform coverage
    data['log_coverage'] = np.log10(data['coverage'].replace(0, np.nan))
    
    # Plot scatter with size based on contig length
    for phylum in unique_phyla:
        subset = data[data['phylum'] == phylum]
        if len(subset) > 0:
            # Scale sizes
            sizes = (subset['length'] / data['length'].max()) * 2000 + 10
            ax_blob.scatter(subset['GC'], subset['log_coverage'],
                          s=sizes, c=colors[phylum], alpha=0.8,
                          edgecolors='gray', linewidths=0.5, label=phylum)
    
    ax_blob.set_xlabel('GC content', fontsize=25)
    ax_blob.set_ylabel('Coverage (log)', fontsize=25)
    ax_blob.set_xlim(0, 1)
    ax_blob.set_ylim(0, 4)
    ax_blob.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x*100)}%'))
    
    # Set y-axis to show powers of 10
    ax_blob.set_yticks([0, 1, 2, 3, 4])
    ax_blob.set_yticklabels(['', '$10^1$', '$10^2$', '$10^3$', '$10^4$'])
    
    ax_blob.tick_params(axis='both', labelsize=20)
    ax_blob.grid(True, alpha=0.3)
    ax_blob.legend(title='Phylum', fontsize=20, title_fontsize=22,
                  loc='upper right', frameon=True)
    
    # Bottom panel: Bar plot
    ax_bar = fig.add_subplot(gs[1])
    
    x_pos = np.arange(len(phylum_counts))
    bars = ax_bar.bar(x_pos, phylum_counts['percent'],
                      color=[colors[p] for p in phylum_counts['phylum']],
                      edgecolor='black', linewidth=1)
    
    # Add percentage labels
    for i, row in phylum_counts.iterrows():
        ax_bar.text(i, row['percent'], f"{row['percent']}%",
                   ha='center', va='bottom', fontsize=20)
    
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(phylum_counts['phylum'], rotation=40, ha='right', fontsize=20)
    ax_bar.set_ylim(0, 100)
    ax_bar.set_ylabel('Percentage', fontsize=25)
    ax_bar.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y)}%'))
    ax_bar.tick_params(axis='y', labelsize=18)
    ax_bar.grid(axis='y', linestyle='--', alpha=0.7)
    ax_bar.set_axisbelow(True)
    
    plt.savefig('Phylum_blob_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save figure object as pickle (Python equivalent of RDS)
    with open('Phylum_final_blobplot.pkl', 'wb') as f:
        pickle.dump(fig, f)
    
    print("Successfully generated Phylum_blob_plot.png")
    print(f"  - Total contigs: {len(data)}")
    print(f"  - Phyla: {len(unique_phyla)}")
    print(f"  - Phylum distribution: {dict(zip(phylum_counts['phylum'], phylum_counts['percent']))}")


if __name__ == '__main__':
    main()
