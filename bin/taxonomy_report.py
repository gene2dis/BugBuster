#!/usr/bin/env python3
"""
Unified Taxonomy Report Generator for BugBuster Pipeline

Generates classification reports and plots for both Kraken2 and Sourmash taxonomic profilers.
Replaces Tax_unify_report.R and SM_unify_report.R with improved performance and maintainability.

Author: BugBuster Development Team
License: MIT
"""

import argparse
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


class TaxonomyReportGenerator:
    """Generate unified taxonomy classification reports and visualizations."""
    
    def __init__(self, profiler: str, db_name: str, output_dir: Path):
        """
        Initialize the report generator.
        
        Args:
            profiler: Taxonomic profiler used ('kraken2' or 'sourmash')
            db_name: Database name for labeling
            output_dir: Output directory for reports and plots
        """
        self.profiler = profiler.lower()
        self.db_name = db_name
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if self.profiler not in ['kraken2', 'sourmash']:
            raise ValueError(f"Profiler must be 'kraken2' or 'sourmash', got '{profiler}'")
    
    def parse_kraken2_reports(self, report_files: List[Path]) -> pd.DataFrame:
        """
        Parse Kraken2 report TSV files.
        
        Expected format: Id, Kraken DB, Unclassified, Classified
        
        Args:
            report_files: List of Kraken2 report TSV files
            
        Returns:
            DataFrame with combined Kraken2 reports
        """
        all_reports = []
        
        for report_file in report_files:
            try:
                df = pd.read_csv(
                    report_file,
                    sep='\t',
                    dtype={
                        'Id': str,
                        'Kraken DB': str,
                        'Unclassified': float,
                        'Classified': float
                    }
                )
                all_reports.append(df)
            except Exception as e:
                print(f"Warning: Failed to parse {report_file}: {e}", file=sys.stderr)
                continue
        
        if not all_reports:
            raise ValueError("No valid Kraken2 reports found")
        
        combined = pd.concat(all_reports, ignore_index=True)
        return combined
    
    def parse_sourmash_reports(self, report_files: List[Path]) -> pd.DataFrame:
        """
        Parse Sourmash report TSV files.
        
        Expected format: Id, Sourmash DB, Unclassified, Classified
        
        Args:
            report_files: List of Sourmash report TSV files
            
        Returns:
            DataFrame with combined Sourmash reports
        """
        all_reports = []
        
        for report_file in report_files:
            try:
                df = pd.read_csv(
                    report_file,
                    sep='\t',
                    dtype={
                        'Id': str,
                        'Sourmash DB': str,
                        'Unclassified': float,
                        'Classified': float
                    }
                )
                all_reports.append(df)
            except Exception as e:
                print(f"Warning: Failed to parse {report_file}: {e}", file=sys.stderr)
                continue
        
        if not all_reports:
            raise ValueError("No valid Sourmash reports found")
        
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
        except FileNotFoundError:
            print(f"Warning: Reads report not found at {reads_report_path}", file=sys.stderr)
            return taxonomy_data
        except Exception as e:
            print(f"Warning: Failed to read reads report: {e}", file=sys.stderr)
            return taxonomy_data
        
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
        if self.profiler == 'kraken2':
            db_col = 'Kraken DB'
            plot_title_prefix = 'Kraken2'
        else:
            db_col = 'Sourmash DB'
            plot_title_prefix = 'Sourmash'
        
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
        
        # Create figure
        fig, ax = plt.subplots(figsize=(16, 8))
        
        # Get sample IDs and sort
        samples = plot_data['Id'].unique()
        y_positions = range(len(samples))
        
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
        ax.set_yticklabels(samples, fontsize=18)
        ax.set_xlabel('Reads Classified', fontsize=20)
        ax.set_ylabel('Samples', fontsize=20)
        
        # Format x-axis as percentage
        ax.set_xlim(0, 1.0)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=18)
        
        # Add legend
        ax.legend(loc='upper right', fontsize=16, frameon=True)
        
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
            taxonomy_data = self.parse_sourmash_reports(report_files)
        
        print(f"Parsed {len(taxonomy_data)} records from {len(report_files)} report files")
        
        # Merge with reads report
        merged_data = self.merge_with_reads_report(taxonomy_data, reads_report_path)
        
        # Save updated reads report
        output_csv = self.output_dir / 'Reads_report.csv'
        merged_data.to_csv(output_csv, index=False)
        print(f"Saved merged report: {output_csv}")
        
        # Create classification plot
        if self.profiler == 'kraken2':
            plot_name = 'Kraken_plot.png'
        else:
            plot_name = 'sourmash_tax_classified_reads.png'
        
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
        """
    )
    
    parser.add_argument(
        '--profiler',
        required=True,
        choices=['kraken2', 'sourmash'],
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
