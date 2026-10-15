#!/usr/bin/env python3
"""
Unified Taxonomy Report Generator for BugBuster Pipeline

Generates classification reports and plots for the Kraken2, Sourmash and Sylph taxonomic profilers.
Replaces Tax_unify_report.R and SM_unify_report.R with improved performance and maintainability.

Author: BugBuster Development Team
License: MIT
"""

import argparse
import re
import sys
from pathlib import Path
from typing import List, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams

# Set matplotlib parameters for publication-quality plots
rcParams['font.size'] = 20
rcParams['axes.linewidth'] = 1.5
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.width'] = 1.5

# Per profiler: plot title label, database column in Reads_report.csv, plot file
PROFILERS = {
    'kraken2': ('Kraken2', 'Kraken DB', 'Kraken_plot.png'),
    'sourmash': ('Sourmash', 'Sourmash DB', 'sourmash_tax_classified_reads.png'),
    'sylph': ('Sylph', 'Sylph DB', 'sylph_tax_classified_reads.png'),
}


def natural_key(value) -> list:
    """Sort key that orders V2 before V10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(value))]


class TaxonomyReportGenerator:
    """Generate unified taxonomy classification reports and visualizations."""
    
    def __init__(self, profiler: str, db_name: str, output_dir: Path):
        """
        Initialize the report generator.
        
        Args:
            profiler: Taxonomic profiler used ('kraken2', 'sourmash' or 'sylph')
            db_name: Database name for labeling
            output_dir: Output directory for reports and plots
        """
        self.profiler = profiler.lower()
        self.db_name = db_name
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if self.profiler not in PROFILERS:
            raise ValueError(f"Profiler must be one of {', '.join(PROFILERS)}, got '{profiler}'")
    
    def parse_kraken2_reports(self, report_files: List[Path]) -> pd.DataFrame:
        """
        Parse standard kraken2 report files (nf-core KRAKEN2 output).

        Each report is a headerless 6-column TSV: percent of reads, clade
        reads, taxon reads, rank code, taxid, name. The sample id comes from
        the filename (``<id>.kraken2.report.txt``); the unclassified fraction
        is the percent of the ``U`` row (0 when every read was classified).

        Args:
            report_files: List of raw kraken2 report files

        Returns:
            DataFrame with columns Id, Kraken DB, Unclassified, Classified
        """
        rows = []

        for report_file in report_files:
            try:
                df = pd.read_csv(
                    report_file,
                    sep='\t',
                    header=None,
                    names=['percent', 'clade_reads', 'taxon_reads', 'rank', 'taxid', 'name'],
                )
            except Exception as e:
                raise ValueError(f"Failed to parse kraken2 report {report_file}: {e}")

            if df.empty or df[['rank', 'taxid', 'name']].isna().any().any():
                raise ValueError(
                    f"{report_file} does not look like a kraken2 report "
                    "(expected 6 tab-separated columns: percent, clade reads, "
                    "taxon reads, rank code, taxid, name)"
                )

            percent = pd.to_numeric(df['percent'], errors='coerce')
            if percent.isna().any():
                raise ValueError(f"Non-numeric percent column in kraken2 report {report_file}")

            unclassified_rows = percent[df['rank'] == 'U']
            unclassified = float(unclassified_rows.iloc[0]) if len(unclassified_rows) else 0.0

            sample_id = re.sub(r'(\.kraken2)?(\.report)?\.txt$', '', report_file.name)

            rows.append({
                'Id': sample_id,
                'Kraken DB': self.db_name,
                'Unclassified': unclassified,
                'Classified': 100.0 - unclassified,
            })

        if not rows:
            raise ValueError("No valid Kraken2 reports found")

        return pd.DataFrame(rows, columns=['Id', 'Kraken DB', 'Unclassified', 'Classified'])
    
    def parse_per_sample_reports(self, report_files: List[Path]) -> pd.DataFrame:
        """
        Parse the per-sample report TSV files written by the SOURMASH and
        SYLPH_PROFILE modules.
        
        Expected format: Id, <Tool> DB, Unclassified, Classified (fractions);
        any further column (SYLPH_PROFILE's 'Read pairs', used by
        taxonomy_phyloseq.py) is dropped so it does not reach Reads_report.csv.
        A report missing one of the four columns raises: skipping it would
        leave that sample's classification blank in Reads_report.csv
        
        Args:
            report_files: List of per-sample report TSV files
            
        Returns:
            DataFrame with combined reports
        """
        label, db_col, _plot_name = PROFILERS[self.profiler]
        columns = ['Id', db_col, 'Unclassified', 'Classified']
        all_reports = []
        
        for report_file in report_files:
            try:
                df = pd.read_csv(
                    report_file,
                    sep='\t',
                    dtype={
                        'Id': str,
                        db_col: str,
                        'Unclassified': float,
                        'Classified': float
                    }
                )
            except Exception as e:
                print(f"Warning: Failed to parse {report_file}: {e}", file=sys.stderr)
                continue
            missing = [c for c in columns if c not in df.columns]
            if missing:
                raise ValueError(f"{report_file} is missing report columns {missing} (expected {columns})")
            all_reports.append(df[columns])
        
        if not all_reports:
            raise ValueError(f"No valid {label} reports found")
        
        combined = pd.concat(all_reports, ignore_index=True)
        return combined
    
    def merge_with_reads_report(
        self,
        taxonomy_data: pd.DataFrame,
        reads_report_path: Path
    ) -> pd.DataFrame:
        """
        Merge taxonomy data with existing reads report.
        
        Args:
            taxonomy_data: DataFrame with taxonomy classification data
            reads_report_path: Path to existing Reads_report.csv
            
        Returns:
            Merged DataFrame
        """
        try:
            reads_report = pd.read_csv(reads_report_path)
        except Exception as e:
            raise ValueError(
                f"Required reads report {reads_report_path} is missing or "
                f"unreadable: {e}"
            )
        
        # Merge on 'Id' column
        # Use suffixes to handle any overlapping columns
        merged = pd.merge(
            reads_report,
            taxonomy_data,
            on='Id',
            how='left',
            suffixes=('', '_taxonomy')
        )
        
        # Drop duplicate columns that came from taxonomy_data (keep the ones from reads_report)
        cols_to_drop = [col for col in merged.columns if col.endswith('_taxonomy')]
        if cols_to_drop:
            merged = merged.drop(columns=cols_to_drop)
        
        return merged
    
    def create_classification_plot(
        self,
        taxonomy_data: pd.DataFrame,
        output_path: Path
    ) -> None:
        """
        Create stacked bar plot showing classified/unclassified reads.
        
        Args:
            taxonomy_data: DataFrame with classification data
            output_path: Path for output PNG file
        """
        # Prepare data for plotting
        plot_title_prefix, db_col, _plot_name = PROFILERS[self.profiler]
        
        # Convert to long format for plotting
        plot_data = taxonomy_data.copy()
        
        # Ensure we have the required columns
        if 'Unclassified' not in plot_data.columns or 'Classified' not in plot_data.columns:
            print(f"Warning: Missing classification columns, skipping plot", file=sys.stderr)
            return
        
        # Convert percentages to proper format (0-1 scale if needed)
        if plot_data['Classified'].max() > 1.0:
            plot_data['Classified'] = plot_data['Classified'] / 100.0
            plot_data['Unclassified'] = plot_data['Unclassified'] / 100.0
        
        # Samples top to bottom in natural order (V2 before V10); figure height
        # and label size follow the number of samples so labels never overlap
        samples = sorted(plot_data['Id'].unique(), key=natural_key)
        n_samples = len(samples)
        fig, ax = plt.subplots(figsize=(16, max(6, 0.3 * n_samples + 2)))
        y_positions = range(n_samples)
        
        # Plot stacked bars
        unclassified_vals = []
        classified_vals = []
        
        for sample in samples:
            sample_data = plot_data[plot_data['Id'] == sample]
            if len(sample_data) > 0:
                unclassified_vals.append(sample_data['Unclassified'].iloc[0])
                classified_vals.append(sample_data['Classified'].iloc[0])
            else:
                unclassified_vals.append(0)
                classified_vals.append(0)
        
        # Create stacked horizontal bars
        ax.barh(y_positions, classified_vals, color='#6baed6', 
                edgecolor='black', linewidth=1.5, label='Classified')
        ax.barh(y_positions, unclassified_vals, left=classified_vals, 
                color='#d16678', edgecolor='black', linewidth=1.5, label='Unclassified')
        
        # Customize plot
        ax.set_yticks(y_positions)
        ax.set_yticklabels(samples, fontsize=max(6, min(18, 500 / max(n_samples, 1))))
        ax.set_ylim(n_samples - 0.5, -0.5)
        ax.set_xlabel('Reads Classified', fontsize=20)
        ax.set_ylabel('Samples', fontsize=20)
        
        # Format x-axis as percentage
        ax.set_xlim(0, 1.0)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=18)
        
        # Add legend
        # outside the axes: inside, it covered the bars of the top samples
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=16, frameon=False)
        
        # Add grid
        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Add title if multiple databases
        if db_col in plot_data.columns:
            db_names = plot_data[db_col].unique()
            if len(db_names) > 1:
                title = f'{plot_title_prefix} Classification Results'
            else:
                title = f'{plot_title_prefix} Classification Results - {db_names[0].upper()}'
            ax.set_title(title, fontsize=22, pad=20)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Created classification plot: {output_path}")
    
    def generate_report(
        self,
        report_files: List[Path],
        reads_report_path: Path
    ) -> Tuple[pd.DataFrame, Path]:
        """
        Generate complete taxonomy report with plots.
        
        Args:
            report_files: List of taxonomy report TSV files
            reads_report_path: Path to existing reads report
            
        Returns:
            Tuple of (merged_dataframe, plot_path)
        """
        # Parse reports based on profiler type
        if self.profiler == 'kraken2':
            taxonomy_data = self.parse_kraken2_reports(report_files)
        else:
            taxonomy_data = self.parse_per_sample_reports(report_files)
        
        print(f"Parsed {len(taxonomy_data)} records from {len(report_files)} report files")
        
        # Merge with reads report
        merged_data = self.merge_with_reads_report(taxonomy_data, reads_report_path)
        
        # Save updated reads report
        output_csv = self.output_dir / 'Reads_report.csv'
        merged_data.to_csv(output_csv, index=False)
        print(f"Saved merged report: {output_csv}")
        
        # Create classification plot
        plot_name = PROFILERS[self.profiler][2]
        
        plot_path = self.output_dir / plot_name
        
        try:
            self.create_classification_plot(taxonomy_data, plot_path)
        except Exception as e:
            print(f"Warning: Failed to create plot: {e}", file=sys.stderr)
            plot_path = None
        
        return merged_data, plot_path


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate unified taxonomy classification reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Kraken2 reports
  taxonomy_report.py --profiler kraken2 --reports *.report.tsv \\
      --reads-report Reads_report.csv --db-name silva --output-dir .

  # Sourmash reports
  taxonomy_report.py --profiler sourmash --reports *.report.tsv \\
      --reads-report Reads_report.csv --db-name gtdb --output-dir .

  # Sylph reports (SYLPH_PROFILE *_report.tsv)
  taxonomy_report.py --profiler sylph --reports *_report.tsv \\
      --reads-report Reads_report.csv --db-name gtdb_r232 --output-dir .
        """
    )
    
    parser.add_argument(
        '--profiler',
        required=True,
        choices=list(PROFILERS),
        help='Taxonomic profiler used'
    )
    
    parser.add_argument(
        '--reports',
        required=True,
        nargs='+',
        type=Path,
        help='Taxonomy report TSV files'
    )
    
    parser.add_argument(
        '--reads-report',
        required=True,
        type=Path,
        help='Path to existing Reads_report.csv'
    )
    
    parser.add_argument(
        '--db-name',
        required=True,
        help='Database name for labeling'
    )
    
    parser.add_argument(
        '--output-dir',
        required=True,
        type=Path,
        help='Output directory for reports and plots'
    )
    
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Validate input files
    valid_reports = [f for f in args.reports if f.exists()]
    if not valid_reports:
        print("Error: No valid report files found", file=sys.stderr)
        sys.exit(1)
    
    if len(valid_reports) < len(args.reports):
        missing = len(args.reports) - len(valid_reports)
        print(f"Warning: {missing} report file(s) not found", file=sys.stderr)
    
    # Create report generator
    generator = TaxonomyReportGenerator(
        profiler=args.profiler,
        db_name=args.db_name,
        output_dir=args.output_dir
    )
    
    # Generate report
    try:
        merged_data, plot_path = generator.generate_report(
            report_files=valid_reports,
            reads_report_path=args.reads_report
        )
        
        print(f"\nReport generation complete!")
        print(f"  - Processed {len(merged_data)} samples")
        print(f"  - Output: {args.output_dir / 'Reads_report.csv'}")
        if plot_path:
            print(f"  - Plot: {plot_path}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
