#!/usr/bin/env python3
"""
Phyloseq-Compatible Table Generator for BugBuster Pipeline

Generates phyloseq-compatible tables and visualizations from Kraken2/Bracken and Sourmash outputs.
Replaces Tax_kraken_to_phyloseq.R, Tax_sourmash_to_phyloseq.R, and eliminates KRAKEN_BIOM process.

Outputs:
- OTU/abundance tables (TSV)
- Taxonomy tables (TSV)
- Sample metadata (TSV)
- HDF5 format for Python analysis
- Abundance plots at multiple taxonomic levels

Author: BugBuster Development Team
License: MIT
"""

import argparse
import sys
import json
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set matplotlib parameters
rcParams['font.size'] = 20
rcParams['axes.linewidth'] = 1.5
rcParams['figure.dpi'] = 300


class PhyloseqTableGenerator:
    """Generate phyloseq-compatible tables from taxonomy profiling data."""
    
    def __init__(
        self,
        profiler: str,
        db_name: str,
        output_dir: Path,
        plot_levels: List[str],
        top_n: int = 10
    ):
        """
        Initialize the phyloseq table generator.
        
        Args:
            profiler: Taxonomic profiler ('kraken2' or 'sourmash')
            db_name: Database name for labeling
            output_dir: Output directory
            plot_levels: Taxonomic levels to plot
            top_n: Number of top taxa to show in plots
        """
        self.profiler = profiler.lower()
        self.db_name = db_name
        self.output_dir = Path(output_dir)
        self.plot_levels = plot_levels
        self.top_n = top_n
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plot_dir = self.output_dir / 'plots'
        self.plot_dir.mkdir(exist_ok=True)
        
        if self.profiler not in ['kraken2', 'sourmash']:
            raise ValueError(f"Profiler must be 'kraken2' or 'sourmash', got '{profiler}'")
        
        self.tax_ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    def parse_biom_file(self, biom_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parse BIOM format file (Kraken2/Bracken output).
        
        Args:
            biom_path: Path to BIOM file
            
        Returns:
            Tuple of (otu_table, taxonomy_table)
        """
        try:
            with open(biom_path, 'r') as f:
                biom_data = json.load(f)
        except Exception as e:
            raise ValueError(f"Failed to parse BIOM file {biom_path}: {e}")
        
        # Extract matrix data
        matrix_type = biom_data.get('matrix_type', 'sparse')
        shape = biom_data.get('shape', [0, 0])
        
        # Initialize abundance matrix
        n_taxa = shape[0]
        n_samples = shape[1]
        abundance_matrix = np.zeros((n_taxa, n_samples))
        
        # Fill matrix based on type
        if matrix_type == 'sparse':
            for entry in biom_data.get('data', []):
                row_idx, col_idx, value = entry
                abundance_matrix[row_idx, col_idx] = value
        else:
            abundance_matrix = np.array(biom_data.get('data', []))
        
        # Extract row (taxa) information
        taxa_ids = []
        taxa_metadata = []
        
        for row in biom_data.get('rows', []):
            taxa_ids.append(row['id'])
            metadata = row.get('metadata', {})
            taxonomy = metadata.get('taxonomy', [])
            taxa_metadata.append(taxonomy)
        
        # Extract column (sample) information
        sample_ids = []
        for col in biom_data.get('columns', []):
            sample_ids.append(col['id'])
        
        # Create OTU table
        otu_table = pd.DataFrame(
            abundance_matrix,
            index=taxa_ids,
            columns=sample_ids
        )
        
        # Create taxonomy table
        tax_table = self._parse_taxonomy_list(taxa_ids, taxa_metadata)
        
        return otu_table, tax_table
    
    def _parse_taxonomy_list(
        self,
        taxa_ids: List[str],
        taxa_metadata: List[List[str]]
    ) -> pd.DataFrame:
        """
        Parse taxonomy metadata into structured table.
        
        Args:
            taxa_ids: List of taxon IDs
            taxa_metadata: List of taxonomy lineages
            
        Returns:
            DataFrame with taxonomic ranks as columns
        """
        tax_data = []
        
        for taxon_id, lineage in zip(taxa_ids, taxa_metadata):
            tax_dict = {'taxon_id': taxon_id}
            
            # Parse lineage (format: "d__Bacteria; p__Proteobacteria; ...")
            for i, rank in enumerate(self.tax_ranks):
                if i < len(lineage):
                    # Remove rank prefix (e.g., "d__", "p__", etc.)
                    tax_value = lineage[i]
                    if '__' in tax_value:
                        tax_value = tax_value.split('__', 1)[1]
                    tax_dict[rank] = tax_value
                else:
                    tax_dict[rank] = ''
            
            tax_data.append(tax_dict)
        
        tax_table = pd.DataFrame(tax_data)
        tax_table.set_index('taxon_id', inplace=True)
        
        return tax_table
    
    def parse_sourmash_gather(
        self,
        gather_files: List[Path]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parse Sourmash gather CSV files with lineage information.
        
        Args:
            gather_files: List of Sourmash gather CSV files
            
        Returns:
            Tuple of (otu_table, taxonomy_table)
        """
        all_data = []
        
        for gather_file in gather_files:
            try:
                df = pd.read_csv(gather_file)
                
                # Extract sample name from filename
                sample_name = gather_file.stem.replace('_smgather', '').replace(f'_{self.db_name}', '')
                df['sample_id'] = sample_name
                
                all_data.append(df)
            except Exception as e:
                print(f"Warning: Failed to parse {gather_file}: {e}", file=sys.stderr)
                continue
        
        if not all_data:
            raise ValueError("No valid Sourmash gather files found")
        
        combined = pd.concat(all_data, ignore_index=True)
        
        # Calculate abundance (unique k-mers * average abundance)
        combined['n_unique_kmers'] = (
            combined['unique_intersect_bp'] / combined['scaled']
        ) * combined['average_abund']
        
        # Extract genome name (remove strain info)
        combined['name'] = combined['name'].str.replace(r' .*', '', regex=True)
        
        # Create OTU table (taxa × samples)
        otu_table = combined.pivot_table(
            index='name',
            columns='sample_id',
            values='n_unique_kmers',
            fill_value=0
        )
        
        # Create taxonomy table
        tax_data = []
        for name, lineage in combined[['name', 'lineage']].drop_duplicates().values:
            tax_dict = {'taxon_id': name}
            
            # Parse lineage (format: "d__Bacteria;p__Proteobacteria;...")
            lineage_parts = lineage.split(';')
            
            for i, rank in enumerate(self.tax_ranks):
                if i < len(lineage_parts):
                    tax_value = lineage_parts[i]
                    # Remove rank prefix
                    if '__' in tax_value:
                        tax_value = tax_value.split('__', 1)[1]
                    tax_dict[rank] = tax_value
                else:
                    tax_dict[rank] = ''
            
            tax_data.append(tax_dict)
        
        tax_table = pd.DataFrame(tax_data)
        tax_table.set_index('taxon_id', inplace=True)
        
        return otu_table, tax_table
    
    def create_sample_metadata(
        self,
        otu_table: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Create sample metadata table.
        
        Args:
            otu_table: OTU abundance table
            
        Returns:
            DataFrame with sample metadata
        """
        metadata = []
        
        for sample_id in otu_table.columns:
            total_reads = otu_table[sample_id].sum()
            n_taxa = (otu_table[sample_id] > 0).sum()
            
            metadata.append({
                'sample_id': sample_id,
                'profiler': self.profiler,
                'database': self.db_name,
                'total_abundance': total_reads,
                'n_taxa_detected': n_taxa
            })
        
        metadata_df = pd.DataFrame(metadata)
        metadata_df.set_index('sample_id', inplace=True)
        
        return metadata_df
    
    def save_tables(
        self,
        otu_table: pd.DataFrame,
        tax_table: pd.DataFrame,
        sample_metadata: pd.DataFrame,
        save_hdf5: bool = True
    ) -> Dict[str, Path]:
        """
        Save phyloseq-compatible tables.
        
        Args:
            otu_table: OTU abundance table
            tax_table: Taxonomy table
            sample_metadata: Sample metadata
            save_hdf5: Whether to save HDF5 format
            
        Returns:
            Dictionary of output file paths
        """
        output_files = {}
        
        # Save TSV tables
        prefix = f"{self.profiler}_{self.db_name}"
        
        otu_path = self.output_dir / f"{prefix}_otu_table.tsv"
        otu_table.to_csv(otu_path, sep='\t')
        output_files['otu_table'] = otu_path
        print(f"Saved OTU table: {otu_path}")
        
        tax_path = self.output_dir / f"{prefix}_tax_table.tsv"
        tax_table.to_csv(tax_path, sep='\t')
        output_files['tax_table'] = tax_path
        print(f"Saved taxonomy table: {tax_path}")
        
        metadata_path = self.output_dir / f"{prefix}_sample_metadata.tsv"
        sample_metadata.to_csv(metadata_path, sep='\t')
        output_files['sample_metadata'] = metadata_path
        print(f"Saved sample metadata: {metadata_path}")
        
        # Save HDF5 format
        if save_hdf5:
            h5_path = self.output_dir / f"{prefix}_phyloseq_data.h5"
            with h5py.File(h5_path, 'w') as f:
                # Store OTU table
                f.create_dataset('otu_table', data=otu_table.values)
                f.create_dataset('otu_taxa', data=[str(x).encode('utf-8') for x in otu_table.index])
                f.create_dataset('otu_samples', data=[str(x).encode('utf-8') for x in otu_table.columns])
                
                # Store taxonomy table
                f.create_dataset('tax_table', data=[[str(v).encode('utf-8') for v in row] for row in tax_table.values])
                f.create_dataset('tax_taxa', data=[str(x).encode('utf-8') for x in tax_table.index])
                f.create_dataset('tax_ranks', data=[str(x).encode('utf-8') for x in tax_table.columns])
                
                # Store metadata
                f.create_dataset('metadata', data=[[str(v).encode('utf-8') for v in row] for row in sample_metadata.values])
                f.create_dataset('metadata_samples', data=[str(x).encode('utf-8') for x in sample_metadata.index])
                f.create_dataset('metadata_columns', data=[str(x).encode('utf-8') for x in sample_metadata.columns])
            
            output_files['hdf5'] = h5_path
            print(f"Saved HDF5 format: {h5_path}")
        
        return output_files
    
    def create_abundance_plot(
        self,
        otu_table: pd.DataFrame,
        tax_table: pd.DataFrame,
        tax_level: str
    ) -> Optional[Path]:
        """
        Create relative abundance bar plot at specified taxonomic level.
        
        Args:
            otu_table: OTU abundance table
            tax_table: Taxonomy table
            tax_level: Taxonomic level to plot
            
        Returns:
            Path to output plot or None if failed
        """
        if tax_level not in tax_table.columns:
            print(f"Warning: Tax level '{tax_level}' not found in taxonomy table", file=sys.stderr)
            return None
        
        try:
            # Aggregate by taxonomic level
            tax_otu = otu_table.copy()
            tax_otu['taxonomy'] = tax_table[tax_level]
            
            # Group by taxonomy
            agg_table = tax_otu.groupby('taxonomy').sum()
            
            # Remove empty/unclassified
            agg_table = agg_table[agg_table.index != '']
            agg_table = agg_table[~agg_table.index.str.contains('unclassified|unknown', case=False, na=False)]
            
            if len(agg_table) == 0:
                print(f"Warning: No classified taxa at {tax_level} level", file=sys.stderr)
                return None
            
            # Calculate relative abundance
            rel_abundance = agg_table.div(agg_table.sum(axis=0), axis=1) * 100
            
            # Select top N taxa
            mean_abundance = rel_abundance.mean(axis=1).sort_values(ascending=False)
            top_taxa = mean_abundance.head(self.top_n).index
            
            plot_data = rel_abundance.loc[top_taxa]
            
            # Add "Other" category
            other = rel_abundance.loc[~rel_abundance.index.isin(top_taxa)].sum(axis=0)
            plot_data.loc['Other'] = other
            
            # Create plot
            fig, ax = plt.subplots(figsize=(16, 8))
            
            # Generate color palette
            n_colors = len(plot_data)
            colors = sns.color_palette("tab20", n_colors)
            
            # Create stacked bar plot
            plot_data.T.plot(
                kind='barh',
                stacked=True,
                ax=ax,
                color=colors,
                edgecolor='black',
                linewidth=0.5,
                width=0.9
            )
            
            # Customize plot
            ax.set_xlabel('Relative Abundance (%)', fontsize=20)
            ax.set_ylabel('Samples', fontsize=20)
            ax.set_title(
                f'Relative Abundance {tax_level} Plot from {self.db_name.upper()}',
                fontsize=22,
                pad=20
            )
            
            # Format x-axis
            ax.set_xlim(0, 100)
            ax.tick_params(axis='both', labelsize=18)
            
            # Legend
            ax.legend(
                bbox_to_anchor=(1.05, 1),
                loc='upper left',
                fontsize=12,
                frameon=True
            )
            
            # Style
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(False)
            
            plt.tight_layout()
            
            # Save plot
            plot_path = self.plot_dir / f"{self.profiler}_{self.db_name}_{tax_level.lower()}_bar.png"
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Created {tax_level} abundance plot: {plot_path}")
            return plot_path
            
        except Exception as e:
            print(f"Warning: Failed to create {tax_level} plot: {e}", file=sys.stderr)
            return None
    
    def generate_all_plots(
        self,
        otu_table: pd.DataFrame,
        tax_table: pd.DataFrame
    ) -> List[Path]:
        """
        Generate abundance plots for all specified taxonomic levels.
        
        Args:
            otu_table: OTU abundance table
            tax_table: Taxonomy table
            
        Returns:
            List of generated plot paths
        """
        plot_paths = []
        
        for tax_level in self.plot_levels:
            plot_path = self.create_abundance_plot(otu_table, tax_table, tax_level)
            if plot_path:
                plot_paths.append(plot_path)
        
        return plot_paths
    
    def process(
        self,
        input_files: List[Path],
        output_format: str = 'both'
    ) -> Dict[str, any]:
        """
        Main processing pipeline.
        
        Args:
            input_files: List of input files (BIOM or Sourmash gather)
            output_format: 'tables', 'hdf5', or 'both'
            
        Returns:
            Dictionary with output information
        """
        print(f"Processing {len(input_files)} input files with {self.profiler} profiler")
        
        # Parse input files based on profiler
        if self.profiler == 'kraken2':
            # For Kraken2, we expect BIOM files (from Bracken)
            if len(input_files) == 1:
                otu_table, tax_table = self.parse_biom_file(input_files[0])
            else:
                # Merge multiple BIOM files
                all_otu = []
                all_tax = []
                for biom_file in input_files:
                    otu, tax = self.parse_biom_file(biom_file)
                    all_otu.append(otu)
                    all_tax.append(tax)
                
                # Combine OTU tables
                otu_table = pd.concat(all_otu, axis=1)
                # Use first taxonomy table (should be consistent)
                tax_table = all_tax[0]
        else:
            # For Sourmash, parse gather CSV files
            otu_table, tax_table = self.parse_sourmash_gather(input_files)
        
        print(f"Parsed data: {len(otu_table)} taxa × {len(otu_table.columns)} samples")
        
        # Create sample metadata
        sample_metadata = self.create_sample_metadata(otu_table)
        
        # Save tables
        save_hdf5 = output_format in ['hdf5', 'both']
        output_files = self.save_tables(
            otu_table,
            tax_table,
            sample_metadata,
            save_hdf5=save_hdf5
        )
        
        # Generate plots
        plot_paths = self.generate_all_plots(otu_table, tax_table)
        
        return {
            'otu_table': otu_table,
            'tax_table': tax_table,
            'sample_metadata': sample_metadata,
            'output_files': output_files,
            'plots': plot_paths
        }


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate phyloseq-compatible tables from taxonomy profiling data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Kraken2/Bracken BIOM files
  taxonomy_phyloseq.py --profiler kraken2 --input-files *.report \\
      --db-name silva --output-dir . --format both \\
      --plot-levels Phylum,Family,Genus,Species

  # Sourmash gather CSV files
  taxonomy_phyloseq.py --profiler sourmash --input-files *_smgather*.csv \\
      --db-name gtdb --output-dir . --format both \\
      --plot-levels Phylum,Family,Genus
        """
    )
    
    parser.add_argument(
        '--profiler',
        required=True,
        choices=['kraken2', 'sourmash'],
        help='Taxonomic profiler used'
    )
    
    parser.add_argument(
        '--input-files',
        required=True,
        nargs='+',
        type=Path,
        help='Input files (BIOM for kraken2, gather CSV for sourmash)'
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
        help='Output directory'
    )
    
    parser.add_argument(
        '--format',
        default='both',
        choices=['tables', 'hdf5', 'both'],
        help='Output format (default: both)'
    )
    
    parser.add_argument(
        '--plot-levels',
        default='Phylum,Family,Genus,Species',
        help='Comma-separated taxonomic levels to plot (default: Phylum,Family,Genus,Species)'
    )
    
    parser.add_argument(
        '--top-n',
        type=int,
        default=10,
        help='Number of top taxa to show in plots (default: 10)'
    )
    
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Validate input files
    valid_files = [f for f in args.input_files if f.exists()]
    if not valid_files:
        print("Error: No valid input files found", file=sys.stderr)
        sys.exit(1)
    
    if len(valid_files) < len(args.input_files):
        missing = len(args.input_files) - len(valid_files)
        print(f"Warning: {missing} input file(s) not found", file=sys.stderr)
    
    # Parse plot levels
    plot_levels = [level.strip() for level in args.plot_levels.split(',')]
    
    # Create generator
    generator = PhyloseqTableGenerator(
        profiler=args.profiler,
        db_name=args.db_name,
        output_dir=args.output_dir,
        plot_levels=plot_levels,
        top_n=args.top_n
    )
    
    # Process
    try:
        results = generator.process(
            input_files=valid_files,
            output_format=args.format
        )
        
        print(f"\nPhyloseq table generation complete!")
        print(f"  - Taxa: {len(results['otu_table'])}")
        print(f"  - Samples: {len(results['otu_table'].columns)}")
        print(f"  - Output files: {len(results['output_files'])}")
        print(f"  - Plots generated: {len(results['plots'])}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
