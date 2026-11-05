#!/usr/bin/env python3
"""
Bin Quality Report Generator
Creates scatter plots and summary tables for bin quality assessment.
Python migration of Bin_checkm_general_plot.R (R original in git at 4bbde47:bin/r_scripts_temp/, removed in roadmap Phase 14)
"""

import glob
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style as fs  # noqa: E402

fs.setup()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


# Columns of Mag_quality_summary.csv and of the two per-MAG tables
SUMMARY_COLUMNS = ['Binner', 'Quality classify', 'Total MAGs']
TABLE_COLUMNS = ['Sample', 'Binner', 'Bin name', 'Quality classify', 'Completeness', 'Contamination']


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
    
    # No sample has a bin: binning found nothing on this data (SemiBin and
    # Binette finish with empty bin sets). Warn and write header-only tables
    # and a placeholder plot so the run can finish; a missing CheckM2 file
    # above still fails
    if not all_reports:
        print(
            "WARNING: all CheckM2 quality reports are empty (header-only) - no "
            "bins were assessed. Upstream binning produced no bins for any "
            "sample; writing empty quality tables.",
            file=sys.stderr
        )
        pd.DataFrame(columns=SUMMARY_COLUMNS).to_csv('Mag_quality_summary.csv', index=False)
        pd.DataFrame(columns=TABLE_COLUMNS).to_csv('All_Mag_quality_table.csv', index=False)
        pd.DataFrame(columns=TABLE_COLUMNS).to_csv('Refined_Mag_quality_table.csv', index=False)
        fs.placeholder('.', 'Total_bins_quality_plot', 'No bins were assessed')
        print("Generated empty quality tables and a placeholder Total_bins_quality_plot")
        return
    
    # Combine all reports
    all_quality_reports = pd.concat(all_reports, ignore_index=True)
    
    # Classify quality
    all_quality_reports['Quality classify'] = all_quality_reports.apply(classify_quality, axis=1)
    
    # Determine bin type
    all_quality_reports['Bin_type'] = all_quality_reports['Binner'].apply(
        lambda x: 'Refined MAGs' if x == 'Binette' else 'Raw MAGs'
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
    all_data.columns = TABLE_COLUMNS
    assert list(mag_summary.columns) == SUMMARY_COLUMNS
    
    refined_data = all_data[all_data['Binner'] == 'Binette'].copy()
    
    # The plot: CheckM2 completeness against contamination, raw and refined bins side by side
    fig, axes = plt.subplots(1, 2, figsize=fs.figsize('single', 62), sharey=True, layout='constrained')
    colors = {'Low': '#9A9A9A', 'Mid': fs.CATEGORY[1], 'High': fs.CATEGORY[2]}
    # a few very contaminated bins would squeeze every other bin against the axis: above 25 % they are drawn
    # at the edge as triangles
    worst = float(quality_data['Contamination'].max())
    xmax = 25.0 if worst > 25 else max(10.0, worst * 1.05)
    for ax, bin_type in zip(axes, ['Raw MAGs', 'Refined MAGs']):
        data_subset = quality_data[quality_data['Bin_type'] == bin_type]
        if len(data_subset) == 0:
            ax.text(0.5, 0.5, f'No {bin_type.lower()}', ha='center', va='center', transform=ax.transAxes)
        for quality in ['Low', 'Mid', 'High']:
            subset = data_subset[data_subset['Quality classify'] == quality]
            if len(subset) > 0:
                # a bin without a contamination value is not drawn (and not taken as contaminated)
                inside = subset['Contamination'] <= xmax
                above = subset['Contamination'] > xmax
                ax.scatter(subset.loc[inside, 'Contamination'], subset.loc[inside, 'Completeness'],
                           color=colors[quality], s=6, lw=0, alpha=0.8)
                ax.scatter([xmax] * int(above.sum()), subset.loc[above, 'Completeness'], marker='>',
                           color=colors[quality], s=10, lw=0, clip_on=False)
        ax.axvline(x=5, color='#555555', linestyle=':', linewidth=0.5, zorder=0)
        for y in (50, 90):
            ax.axhline(y=y, color='#555555', linestyle=':', linewidth=0.5, zorder=0)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, 101)
        ax.set_xlabel('contamination (%)')
        ax.set_title(f'{bin_type.split()[0].lower()} bins ({len(data_subset)})', fontsize=fs.BASE_PT,
                     fontweight='normal')
    axes[0].set_ylabel('completeness (%)')
    legend_labels = ['High (completeness >= 90 %, contamination <= 5 %)',
                     'Mid (completeness 50-90 %, contamination <= 5 %)',
                     'Low (completeness < 50 % or contamination > 5 %)']
    handles = [Line2D([], [], ls='', marker='o', markersize=3, color=colors[q], label=lab)
               for q, lab in zip(['High', 'Mid', 'Low'], legend_labels)]
    if worst > xmax:
        handles.append(Line2D([], [], ls='', marker='>', markersize=3, color='#555555',
                              label='contamination above 25 % (at the edge)'))
    fig.legend(handles=handles, title='Quality class (not MIMAG)', loc='outside lower left', ncol=1,
               alignment='left')
    source = quality_data.rename(columns={'Name': 'Bin name'})[
        ['Sample', 'Binner', 'Bin name', 'Bin_type', 'Quality classify', 'Completeness', 'Contamination']]
    fs.save(fig, '.', 'Total_bins_quality_plot', source)
    print("Successfully generated Total_bins_quality_plot")

    # Write output files
    mag_summary.to_csv('Mag_quality_summary.csv', index=False)
    all_data.to_csv('All_Mag_quality_table.csv', index=False)
    refined_data.to_csv('Refined_Mag_quality_table.csv', index=False)
    
    print("Successfully generated quality reports:")
    print(f"  - Mag_quality_summary.csv ({len(mag_summary)} rows)")
    print(f"  - All_Mag_quality_table.csv ({len(all_data)} MAGs)")
    print(f"  - Refined_Mag_quality_table.csv ({len(refined_data)} refined MAGs)")


if __name__ == '__main__':
    main()
