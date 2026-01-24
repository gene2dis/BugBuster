#!/usr/bin/env python3
"""
ARG Blob Plot Generator
Creates multi-panel blob plots showing contig GC vs coverage colored by ARG class.
Python migration of ARG_blob_plot.R
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import numpy as np
import sys
import pickle


def generate_safe_colors(n):
    """Generate color-blind safe palette similar to rcartocolor Safe."""
    # Safe palette colors (color-blind friendly)
    safe_colors = [
        '#88CCEE', '#CC6677', '#DDCC77', '#117733', '#332288',
        '#AA4499', '#44AA99', '#999933', '#882255', '#661100',
        '#6699CC', '#888888'
    ]
    if n <= len(safe_colors):
        return safe_colors[:n]
    else:
        # If more colors needed, cycle through
        return [safe_colors[i % len(safe_colors)] for i in range(n)]


def main():
    """Main function to generate ARG blob plot."""
    
    # Read contig ARG and taxonomy data
    try:
        contig_arg_data = pd.read_csv('Contig_tax_and_arg_prediction.tsv')
    except FileNotFoundError:
        print("ERROR: Contig_tax_and_arg_prediction.tsv not found", file=sys.stderr)
        sys.exit(1)
    
    # Check if all ARG genes are NA
    all_na = contig_arg_data['ARG_gene'].isna().all()
    
    if all_na:
        # If no ARGs, set everything to "None"
        contig_arg_data['ARG_gene'] = contig_arg_data['ARG_gene'].fillna('None')
        contig_arg_data['predicted_ARG_class'] = 'None'
    
    # Calculate total length by ARG class
    arg_class_summary = contig_arg_data.groupby('predicted_ARG_class')['contig_length'].sum().reset_index()
    arg_class_summary.columns = ['predicted_ARG_class', 'length_sum']
    arg_class_summary = arg_class_summary.sort_values('length_sum', ascending=False)
    
    # Determine max number of classes to show (top 9 or all if less)
    max_num_pull = min(9, len(arg_class_summary)) if not all_na else 1
    
    # Prepare data for plotting
    if all_na:
        data_to_plot = contig_arg_data[['contig_name', 'contig_length', 'GC', 'coverage', 'predicted_ARG_class']].copy()
        data_to_plot['predicted_ARG_class'] = data_to_plot['predicted_ARG_class'].str.title()
        data_to_plot['coverage'] = np.log10(data_to_plot['coverage'].replace(0, np.nan))
    else:
        # Keep top classes, group rest as "Other"
        top_classes = arg_class_summary['predicted_ARG_class'].head(max_num_pull).tolist()
        data_to_plot = contig_arg_data[['contig_name', 'contig_length', 'GC', 'coverage', 'predicted_ARG_class']].copy()
        data_to_plot['predicted_ARG_class'] = data_to_plot['predicted_ARG_class'].apply(
            lambda x: x if x in top_classes else 'Other'
        )
        data_to_plot['predicted_ARG_class'] = data_to_plot['predicted_ARG_class'].str.title()
        data_to_plot = data_to_plot.dropna(subset=['predicted_ARG_class'])
        
        # Log transform coverage, handle -Inf
        data_to_plot['coverage'] = np.log10(data_to_plot['coverage'].replace(0, np.nan))
        min_coverage = data_to_plot['coverage'].replace(-np.inf, np.nan).min()
        data_to_plot['coverage'] = data_to_plot['coverage'].replace(-np.inf, min_coverage)
    
    # Generate colors
    unique_classes = sorted(data_to_plot['predicted_ARG_class'].unique())
    if all_na:
        colors = {'None': 'gray'}
    else:
        color_list = generate_safe_colors(len(unique_classes))
        colors = dict(zip(unique_classes, color_list))
    
    # Calculate percentages
    class_counts = data_to_plot.groupby('predicted_ARG_class').size().reset_index(name='n')
    class_counts['percent'] = (class_counts['n'] / class_counts['n'].sum() * 100).round(2)
    
    # Create figure with custom layout
    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.15)
    
    # Top panel: Blob plot
    ax_blob = fig.add_subplot(gs[0])
    
    # Plot scatter with size based on contig length
    for arg_class in unique_classes:
        subset = data_to_plot[data_to_plot['predicted_ARG_class'] == arg_class]
        if len(subset) > 0:
            # Scale sizes
            sizes = (subset['contig_length'] / data_to_plot['contig_length'].max()) * 2000 + 10
            ax_blob.scatter(subset['GC'], subset['coverage'],
                          s=sizes, c=colors[arg_class], alpha=0.8,
                          edgecolors='gray', linewidths=0.5, label=arg_class)
    
    ax_blob.set_xlabel('GC content', fontsize=25)
    ax_blob.set_ylabel('Coverage (log)', fontsize=25)
    ax_blob.set_xlim(0, 1)
    ax_blob.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x*100)}%'))
    ax_blob.tick_params(axis='both', labelsize=20)
    ax_blob.grid(True, alpha=0.3)
    ax_blob.legend(title='Predicted ARG Class', fontsize=20, title_fontsize=22,
                  loc='upper right', frameon=True)
    
    # Bottom panel: Bar plot
    ax_bar = fig.add_subplot(gs[1])
    
    x_pos = np.arange(len(class_counts))
    bars = ax_bar.bar(x_pos, class_counts['percent'],
                      color=[colors[c] for c in class_counts['predicted_ARG_class']],
                      edgecolor='black', linewidth=1)
    
    # Add percentage labels
    for i, row in class_counts.iterrows():
        ax_bar.text(i, row['percent'], f"{row['percent']}%",
                   ha='center', va='bottom', fontsize=20)
    
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(class_counts['predicted_ARG_class'], rotation=40, ha='right', fontsize=20)
    ax_bar.set_ylim(0, 100)
    ax_bar.set_ylabel('Percentage', fontsize=25)
    ax_bar.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y)}%'))
    ax_bar.tick_params(axis='y', labelsize=18)
    ax_bar.grid(axis='y', linestyle='--', alpha=0.7)
    ax_bar.set_axisbelow(True)
    
    plt.savefig('ARG_blob_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save figure object as pickle (Python equivalent of RDS)
    with open('ARG_final_blobplot.pkl', 'wb') as f:
        pickle.dump(fig, f)
    
    print("Successfully generated ARG_blob_plot.png")
    print(f"  - Total contigs: {len(data_to_plot)}")
    print(f"  - ARG classes: {len(unique_classes)}")
    print(f"  - Class distribution: {dict(zip(class_counts['predicted_ARG_class'], class_counts['percent']))}")


if __name__ == '__main__':
    main()
