#!/usr/bin/env python3
"""
Bin Quality Report Generator
Creates scatter plots and summary tables for bin quality assessment.
Python migration of Bin_checkm_general_plot.R
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import glob
import sys
import numpy as np


def classify_quality(row):
    """Classify bin quality based on completeness and contamination."""
    if row['Completeness'] >= 90 and row['Contamination'] <= 5:
        return 'High'
    elif row['Completeness'] >= 50 and row['Contamination'] <= 5:
        return 'Mid'
    else:
        return 'Low'


def main():
    """Main function to generate bin quality report."""
    
    # Find all quality report files
    quality_files = glob.glob('*_quality_report.tsv')
    if not quality_files:
        print("ERROR: No quality report files found (*_quality_report.tsv)", file=sys.stderr)
        sys.exit(1)
    
    # Read and process all quality reports
    all_reports = []
    for file in quality_files:
        # Extract sample name and binner from filename
        # Format: sample_binner_quality_report.tsv
        parts = file.replace('_quality_report.tsv', '').split('_')
        if len(parts) >= 2:
            binner = parts[-1]
            sample = '_'.join(parts[:-1])
        else:
            print(f"WARNING: Could not parse filename {file}, skipping", file=sys.stderr)
            continue
        
        df = pd.read_csv(file, sep='\t')
        
        # Skip empty files (only header, no data rows)
        if len(df) == 0:
            print(f"WARNING: Skipping empty file {file} (no bins assessed)")
            continue
        
        df['Binner'] = binner.capitalize()
        df['Sample'] = sample
        all_reports.append(df)
    
    if not all_reports:
        print("WARNING: No bins were assessed by CheckM2. Creating empty output files.")
        # Create empty output files with proper structure
        empty_summary = pd.DataFrame(columns=['Binner', 'Quality classify', 'Total MAGs'])
        empty_summary.to_csv('Mag_quality_summary.csv', index=False)
        
        empty_table = pd.DataFrame(columns=['Sample', 'Binner', 'Bin name', 'Quality classify', 'Completeness', 'Contamination'])
        empty_table.to_csv('All_Mag_quality_table.csv', index=False)
        empty_table.to_csv('Refined_Mag_quality_table.csv', index=False)
        
        # Create empty plot
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))
        for idx, bin_type in enumerate(['Raw MAGs', 'Refined MAGs']):
            axes[idx].text(0.5, 0.5, f'No {bin_type}', ha='center', va='center', fontsize=20, transform=axes[idx].transAxes)
            axes[idx].set_xlim(0, 20)
            axes[idx].set_ylim(0, 100)
            axes[idx].set_xlabel('Contamination', fontsize=18)
            if idx == 0:
                axes[idx].set_ylabel('Completeness', fontsize=18)
            axes[idx].set_title(bin_type, fontsize=20, pad=10, 
                               bbox=dict(boxstyle='round', facecolor='#CDDEFF', alpha=0.5))
            axes[idx].tick_params(axis='both', labelsize=18)
            axes[idx].yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y)}%'))
            axes[idx].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}%'))
            axes[idx].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('Total_bins_quality_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Created empty output files (no bins to assess)")
        return
    
    # Combine all reports
    all_quality_reports = pd.concat(all_reports, ignore_index=True)
    
    # Classify quality
    all_quality_reports['Quality classify'] = all_quality_reports.apply(classify_quality, axis=1)
    
    # Determine bin type
    all_quality_reports['Bin_type'] = all_quality_reports['Binner'].apply(
        lambda x: 'Refined MAGs' if x == 'Metawrap' else 'Raw MAGs'
    )
    
    # Select relevant columns
    quality_data = all_quality_reports[[
        'Binner', 'Quality classify', 'Completeness', 'Contamination', 
        'Sample', 'Name', 'Bin_type'
    ]].copy()
    
    # Set quality order
    quality_data['Quality classify'] = pd.Categorical(
        quality_data['Quality classify'],
        categories=['Low', 'Mid', 'High'],
        ordered=True
    )
    
    # Generate summary table
    mag_summary = quality_data.groupby(['Binner', 'Quality classify'], observed=False).size().reset_index(name='Total MAGs')
    
    # Prepare data tables for export
    all_data = quality_data[['Sample', 'Binner', 'Name', 'Quality classify', 'Completeness', 'Contamination']].copy()
    all_data.columns = ['Sample', 'Binner', 'Bin name', 'Quality classify', 'Completeness', 'Contamination']
    
    refined_data = all_data[all_data['Binner'] == 'Metawrap'].copy()
    
    # Create the plot
    try:
        fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)
        
        # Define colors
        colors = {'Low': '#d16678', 'Mid': '#6baed6', 'High': '#78b33e'}
        
        # Plot for each bin type
        for idx, bin_type in enumerate(['Raw MAGs', 'Refined MAGs']):
            ax = axes[idx]
            data_subset = quality_data[quality_data['Bin_type'] == bin_type]
            
            if len(data_subset) == 0:
                ax.text(0.5, 0.5, f'No {bin_type}', ha='center', va='center', fontsize=20)
                ax.set_xlim(0, 20)
                ax.set_ylim(0, 100)
            else:
                # Plot points for each quality class
                for quality in ['Low', 'Mid', 'High']:
                    subset = data_subset[data_subset['Quality classify'] == quality]
                    if len(subset) > 0:
                        ax.scatter(subset['Contamination'], subset['Completeness'],
                                 c=colors[quality], s=100, alpha=0.5, label=quality)
                
                # Add reference lines
                ax.axvline(x=5, color='black', linestyle='--', linewidth=0.8)
                ax.axhline(y=90, color='black', linestyle='--', linewidth=0.8)
                ax.axhline(y=50, color='black', linestyle='--', linewidth=0.8)
            
            # Set labels and limits
            ax.set_xlabel('Contamination', fontsize=18)
            if idx == 0:
                ax.set_ylabel('Completeness', fontsize=18)
            ax.set_xlim(0, 20)
            ax.set_ylim(0, 100)
            ax.set_title(bin_type, fontsize=20, pad=10, 
                        bbox=dict(boxstyle='round', facecolor='#CDDEFF', alpha=0.5))
            
            # Format axes
            ax.tick_params(axis='both', labelsize=18)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y)}%'))
            ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}%'))
            ax.grid(True, alpha=0.3)
        
        # Create legend
        legend_labels = [
            'Low\n(Completeness < 50%, Contamination > 5%)',
            'Mid\n(Completeness ≥ 50% < 90%, Contamination ≤ 5%)',
            'High\n(Completeness ≥ 90%, Contamination ≤ 5%)'
        ]
        handles = [plt.Line2D([0], [0], marker='o', color='w', 
                             markerfacecolor=colors[q], markersize=15, alpha=0.5)
                  for q in ['Low', 'Mid', 'High']]
        
        fig.legend(handles, legend_labels, loc='upper right', 
                  fontsize=18, title='Quality classification',
                  title_fontsize=18, frameon=True, bbox_to_anchor=(0.98, 0.98))
        
        plt.tight_layout()
        plt.savefig('Total_bins_quality_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Successfully generated Total_bins_quality_plot.png")
        
    except Exception as e:
        print(f"ERROR generating plot: {e}", file=sys.stderr)
    
    # Write output files
    mag_summary.to_csv('Mag_quality_summary.csv', index=False)
    all_data.to_csv('All_Mag_quality_table.csv', index=False)
    refined_data.to_csv('Refined_Mag_quality_table.csv', index=False)
    
    print(f"Successfully generated quality reports:")
    print(f"  - Mag_quality_summary.csv ({len(mag_summary)} rows)")
    print(f"  - All_Mag_quality_table.csv ({len(all_data)} MAGs)")
    print(f"  - Refined_Mag_quality_table.csv ({len(refined_data)} refined MAGs)")


if __name__ == '__main__':
    main()
