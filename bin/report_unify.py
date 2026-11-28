#!/usr/bin/env python3
"""
Report Unify
Consolidates the QC reports of fastp and the host-removal step (the
decontamination aligner, bwa-mem3 by default or Bowtie2 with --read_aligner
bowtie2; both write the same <id>_<filter>_decontamination_report.tsv).
Generates CSV summary and high-quality scientific boxplot visualization.

Python port of Report_unify.R.
"""

import sys
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style as fs  # noqa: E402
from id_columns import read_table  # noqa: E402

fs.setup()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

warnings.filterwarnings('ignore')


def read_tsv_files(pattern):
    """Read and concatenate TSV files matching pattern."""
    files = sorted(glob.glob(pattern))
    if not files:
        return pd.DataFrame()
    
    dfs = []
    for file in files:
        try:
            # Id as text: one file per sample, so pandas would read 001 as 1
            df = read_table(file, ['Id'], sep='\t')
            dfs.append(df)
        except Exception as e:
            print(f"Warning: Could not read {file}: {e}", file=sys.stderr)
    
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()


def process_reports_with_filtering(filter_args):
    """Process reports when QC filtering was performed."""
    # Read the decontamination reports of each filter type
    decontam_dfs = {}
    for arg in filter_args:
        # the decontamination report's own suffix: a substring match on the
        # whole name also caught a sample named after the alias (R6-09)
        pattern = f"*_{arg}_decontamination_report.tsv"
        df = read_tsv_files(pattern)
        if not df.empty:
            decontam_dfs[f"decontam_{arg}"] = df
    
    # Read fastp reports
    fastp_reports = read_tsv_files("*fastp_report.tsv")
    
    if fastp_reports.empty:
        print("Error: No fastp reports found", file=sys.stderr)
        return None, None
    
    # Merge all decontamination reports
    if decontam_dfs:
        all_decontam = None
        for df in decontam_dfs.values():
            if all_decontam is None:
                all_decontam = df
            else:
                all_decontam = pd.merge(all_decontam, df, on='Id', how='left')
        
        # Merge with fastp reports
        all_reports = pd.merge(fastp_reports, all_decontam, on='Id', how='left')
    else:
        all_reports = fastp_reports
    
    # Fill NaN with 0
    all_reports = all_reports.fillna(0)
    
    # Sort columns by total reads (descending)
    numeric_cols = all_reports.select_dtypes(include=[np.number]).columns
    col_sums = all_reports[numeric_cols].sum().sort_values(ascending=False, kind="stable")
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


def create_boxplot(data_plot, key="Box_plot_reads"):
    """Reads per sample after each step (box and points), one panel per read type, in the shared figure style:
    <key>.pdf / .svg / .png and the plotted values as <key>.csv."""
    if data_plot is None or data_plot.empty:
        print("Warning: No data to plot", file=sys.stderr)
        return
    read_types = sorted(data_plot['Read_type'].unique())
    processes = list(data_plot['Process'].cat.categories)
    pal = dict(zip(processes, (fs.CATEGORY * 3)[:len(processes)]))
    width = "single" if len(read_types) * len(processes) <= 8 else "double"
    fig, axes = plt.subplots(1, len(read_types), figsize=fs.figsize(width, 62), sharey=True, layout="constrained",
                             squeeze=False)
    rng = np.random.default_rng(1)  # fixed jitter: a rerun draws the same figure
    for ax, read_type in zip(axes[0], read_types):
        subset = data_plot[data_plot['Read_type'] == read_type]
        box_data = [subset[subset['Process'] == proc]['Reads'].values for proc in processes]
        bp = ax.boxplot(box_data, positions=range(len(processes)), widths=0.55, patch_artist=True,
                        showfliers=False, boxprops={"linewidth": 0.5}, whiskerprops={"linewidth": 0.5},
                        capprops={"linewidth": 0.5}, medianprops={"linewidth": 0.8, "color": "black"})
        for patch, process in zip(bp['boxes'], processes):
            patch.set_facecolor(pal[process])
            patch.set_alpha(0.35)
        for i, process in enumerate(processes):
            y = subset[subset['Process'] == process]['Reads'].values
            if len(y):
                ax.scatter(i + rng.uniform(-0.12, 0.12, len(y)), y, s=5, lw=0, color=pal[process], zorder=3)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: fs.short_count(y) if y else "0"))
        ax.set_xticks(range(len(processes)), processes, rotation=45, ha="right")
        ax.set_title(str(read_type), fontsize=fs.BASE_PT, fontweight="normal")
    axes[0][0].set_ylabel("reads per sample")
    fs.save(fig, ".", key, data_plot[["Id", "Process", "Read_type", "Reads"]] if "Id" in data_plot else data_plot)
    print(f"Boxplot saved to {key}.png / .pdf / .svg")


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
        if all_reports is None:
            sys.exit(1)
        all_reports.to_csv("Reads_report.csv", index=False)
        print("Report saved to Reads_report.csv (no filtering mode)")
    else:
        # Filtering was performed
        all_reports, data_plot = process_reports_with_filtering(args)
        if all_reports is None:
            sys.exit(1)

        # Save CSV report
        all_reports.to_csv("Reads_report.csv", index=False)
        print("Report saved to Reads_report.csv")

        # Create boxplot
        if data_plot is not None:
            create_boxplot(data_plot)


if __name__ == "__main__":
    main()
