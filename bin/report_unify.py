#!/usr/bin/env python3
"""
Report Unify - Python Migration
Consolidates QC reports from fastp and bowtie2 filtering steps.
Generates CSV summary and high-quality scientific boxplot visualization.

Original R script: Report_unify.R
Migrated to Python for QC subworkflow simplification.
"""

import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def read_tsv_files(pattern):
    """Read and concatenate TSV files matching pattern."""
    files = glob.glob(pattern)
    if not files:
        return pd.DataFrame()
    
    dfs = []
    for file in files:
        try:
            df = pd.read_csv(file, sep='\t')
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {file}: {e}", file=sys.stderr)
    
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()


def process_reports_with_filtering(filter_args):
    """Process reports when QC filtering was performed."""
    # Read bowtie reports for each filter type
    bowtie_dfs = {}
    for arg in filter_args:
        pattern = f"*{arg}*"
        df = read_tsv_files(pattern)
        if not df.empty:
            bowtie_dfs[f"bowtie_{arg}"] = df
    
    # Read fastp reports
    fastp_reports = read_tsv_files("*fastp_report.tsv")
    
    if fastp_reports.empty:
        print("Error: No fastp reports found", file=sys.stderr)
        return None, None
    
    # Merge all bowtie reports
    if bowtie_dfs:
        all_bowtie = None
        for df in bowtie_dfs.values():
            if all_bowtie is None:
                all_bowtie = df
            else:
                all_bowtie = pd.merge(all_bowtie, df, on='Id', how='left')
        
        # Merge with fastp reports
        all_reports = pd.merge(fastp_reports, all_bowtie, on='Id', how='left')
    else:
        all_reports = fastp_reports
    
    # Fill NaN with 0
    all_reports = all_reports.fillna(0)
    
    # Sort columns by total reads (descending)
    numeric_cols = all_reports.select_dtypes(include=[np.number]).columns
    col_sums = all_reports[numeric_cols].sum().sort_values(ascending=False)
    ordered_cols = ['Id'] + col_sums.index.tolist()
    all_reports_order = all_reports[ordered_cols]
    
    # Prepare data for plotting
    data_plot = all_reports_order.melt(
        id_vars=['Id'],
        var_name='Process',
        value_name='Reads'
    )
    
    # Classify read types
    data_plot['Read_type'] = data_plot['Process'].apply(
        lambda x: 'Singleton' if 'singletons' in x.lower() else 'Paired'
    )
    
    # Clean process names
    data_plot['Process'] = data_plot['Process'].str.replace(' singletons', '', regex=False)
    
    # Rename processes for clarity
    def rename_process(process, args):
        if 'Raw' in process:
            return 'Raw reads'
        for i, arg in enumerate(args):
            if arg in process:
                return f"{arg} reads removal"
        return process
    
    data_plot['Process'] = data_plot['Process'].apply(lambda x: rename_process(x, filter_args))
    
    # Set factor levels for ordered plotting
    process_levels = ['Raw reads', 'Fastp'] + [f"{arg} reads removal" for arg in filter_args]
    data_plot['Process'] = pd.Categorical(
        data_plot['Process'],
        categories=process_levels,
        ordered=True
    )
    
    return all_reports_order, data_plot


def process_reports_no_filtering():
    """Process reports when no QC filtering was performed."""
    fastp_reports = read_tsv_files("*fastp_report.tsv")
    
    if fastp_reports.empty:
        print("Error: No fastp reports found", file=sys.stderr)
        return None, None
    
    return fastp_reports, None


def create_boxplot(data_plot, output_file="Box_plot_reads.png"):
    """Create high-quality scientific boxplot matching R ggplot2 style."""
    if data_plot is None or data_plot.empty:
        print("Warning: No data to plot", file=sys.stderr)
        return
    
    try:
        # Calculate y-axis limits
        max_reads = data_plot['Reads'].max()
        max_round = np.ceil(max_reads / 5000000) * 5000000
        
        # Set publication-quality style
        plt.rcParams['font.size'] = 20
        plt.rcParams['axes.linewidth'] = 1.0
        plt.rcParams['xtick.major.size'] = 5
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['xtick.major.width'] = 1.0
        plt.rcParams['ytick.major.width'] = 1.0
        
        # Create figure with facets
        read_types = data_plot['Read_type'].unique()
        n_facets = len(read_types)
        
        fig, axes = plt.subplots(1, n_facets, figsize=(16, 8), sharey=True)
        if n_facets == 1:
            axes = [axes]
        
        # Get unique processes for color palette
        processes = data_plot['Process'].cat.categories
        colors = sns.color_palette("husl", len(processes))
        color_map = dict(zip(processes, colors))
        
        for idx, read_type in enumerate(sorted(read_types)):
            ax = axes[idx]
            subset = data_plot[data_plot['Read_type'] == read_type]
            
            # Create boxplot
            positions = range(len(processes))
            box_data = [subset[subset['Process'] == proc]['Reads'].values 
                       for proc in processes]
            
            bp = ax.boxplot(
                box_data,
                positions=positions,
                widths=0.6,
                patch_artist=True,
                showfliers=True,
                flierprops=dict(marker='o', markerfacecolor='blue', markersize=8, 
                              linestyle='none', markeredgecolor='blue', alpha=0.6),
                boxprops=dict(linewidth=1.5),
                whiskerprops=dict(linewidth=1.5),
                capprops=dict(linewidth=1.5),
                medianprops=dict(linewidth=2, color='black')
            )
            
            # Color boxes
            for patch, process in zip(bp['boxes'], processes):
                patch.set_facecolor(color_map[process])
                patch.set_alpha(0.7)
            
            # Add jittered points
            for i, process in enumerate(processes):
                y_data = subset[subset['Process'] == process]['Reads'].values
                if len(y_data) > 0:
                    x_data = np.random.normal(i, 0.04, size=len(y_data))
                    ax.scatter(x_data, y_data, alpha=0.4, s=80, 
                             color=color_map[process], zorder=3)
            
            # Styling
            ax.set_xlabel('Process', fontsize=20, fontweight='normal')
            if idx == 0:
                ax.set_ylabel('Total reads', fontsize=20, fontweight='normal')
            
            # Y-axis formatting
            y_ticks = np.arange(0, max_round + 5000000, 5000000)
            ax.set_ylim(0, max_round)
            ax.set_yticks(y_ticks)
            ax.set_yticklabels([f'{int(y/1e6)}M' for y in y_ticks])
            
            # X-axis
            ax.set_xticks(positions)
            ax.set_xticklabels(processes, rotation=60, ha='right', va='top')
            
            # Facet label with custom background
            ax.text(0.5, 0.98, read_type, transform=ax.transAxes,
                   fontsize=20, ha='center', va='top',
                   bbox=dict(boxstyle='round,pad=0.5', 
                           facecolor='#CDDEFF', edgecolor='#CDDEFF', linewidth=0))
            
            # Grid and spines
            ax.grid(True, axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
            ax.set_axisbelow(True)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.0)
            ax.spines['bottom'].set_linewidth(1.0)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"Boxplot saved to {output_file}")
        
    except Exception as e:
        print(f"Error creating boxplot: {e}", file=sys.stderr)


def main():
    """Main execution function."""
    if len(sys.argv) < 2:
        print("Usage: report_unify.py <filter_arg1> [filter_arg2] ...", file=sys.stderr)
        print("       report_unify.py none", file=sys.stderr)
        sys.exit(1)
    
    args = sys.argv[1:]
    
    if args[0] == "none":
        # No filtering performed
        all_reports, _ = process_reports_no_filtering()
        if all_reports is not None:
            all_reports.to_csv("Reads_report.csv", index=False)
            print("Report saved to Reads_report.csv (no filtering mode)")
    else:
        # Filtering was performed
        all_reports, data_plot = process_reports_with_filtering(args)
        
        if all_reports is not None:
            # Save CSV report
            all_reports.to_csv("Reads_report.csv", index=False)
            print("Report saved to Reads_report.csv")
            
            # Create boxplot
            if data_plot is not None:
                create_boxplot(data_plot)


if __name__ == "__main__":
    main()
