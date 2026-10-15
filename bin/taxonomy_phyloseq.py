#!/usr/bin/env python3
"""
Phyloseq-Compatible Table Generator for BugBuster Pipeline

Generates phyloseq-compatible tables and visualizations from Kraken2/Bracken, Sourmash and Sylph outputs.
Replaces Tax_kraken_to_phyloseq.R, Tax_sourmash_to_phyloseq.R, and eliminates KRAKEN_BIOM process.

Outputs:
- OTU/abundance tables (TSV)
- Taxonomy tables (TSV)
- Sample metadata (TSV)
- Sylph only: estimated read pairs per genome (TSV, <prefix>_read_counts.tsv)
- HDF5 format for Python analysis
- Abundance plots at multiple taxonomic levels

Author: BugBuster Development Team
License: MIT
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from concurrent.futures import ProcessPoolExecutor, as_completed
import re


def natural_key(value) -> list:
    """Sort key that orders V2 before V10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(value))]

# Set matplotlib parameters
rcParams['font.size'] = 20
rcParams['axes.linewidth'] = 1.5
rcParams['figure.dpi'] = 300

PROFILERS = ['kraken2', 'sourmash', 'sylph']


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
            profiler: Taxonomic profiler ('kraken2', 'sourmash' or 'sylph')
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
        
        if self.profiler not in PROFILERS:
            raise ValueError(f"Profiler must be one of {', '.join(PROFILERS)}, got '{profiler}'")
        
        self.tax_ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    # Rank codes in kraken-style reports mapped to phyloseq ranks. Kraken2
    # emits 'D' (domain) at the top level for NCBI taxonomies while other DB
    # builds use 'K'; both fill the Kingdom slot. Sub-ranks (S1, G1, ...) and
    # unranked rows stay on the indent-walk stack but fill no rank column.
    KRAKEN_RANK_MAP = {
        'D': 'Kingdom', 'K': 'Kingdom', 'P': 'Phylum', 'C': 'Class',
        'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species',
    }

    def parse_kraken_reports(
        self,
        report_files: List[Path]
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parse kraken-style reports (bracken's re-estimated ``-w`` output).

        Each report is a headerless 6-column TSV (percent, clade reads, taxon
        reads, rank code, taxid, name) where the lineage is encoded by rank
        codes plus 2-space name indentation. Species rows become OTUs keyed
        by taxid; their lineage is rebuilt by walking the indentation.

        Args:
            report_files: List of per-sample kraken-style report files

        Returns:
            Tuple of (otu_table, taxonomy_table)
        """
        per_sample_counts = {}
        tax_rows = {}

        for report_file in report_files:
            sample_name = report_file.name
            for suffix in ('.txt', '_bracken', '.report', '.kraken2'):
                if sample_name.endswith(suffix):
                    sample_name = sample_name[: -len(suffix)]

            counts = {}
            stack = []  # (depth, rank_code, name)
            with open(report_file) as fh:
                for line_no, line in enumerate(fh, 1):
                    if not line.strip():
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) != 6:
                        raise ValueError(
                            f"{report_file}:{line_no}: expected 6 tab-separated "
                            f"kraken-report columns, got {len(fields)}"
                        )
                    _percent, clade_reads, _taxon_reads, rank, taxid, name = fields
                    if rank == 'U':
                        continue
                    depth = (len(name) - len(name.lstrip(' '))) // 2
                    while stack and stack[-1][0] >= depth:
                        stack.pop()
                    stack.append((depth, rank, name.strip()))
                    if rank == 'S':
                        try:
                            clade_count = int(clade_reads)
                        except ValueError:
                            raise ValueError(
                                f"{report_file}:{line_no}: non-integer clade "
                                f"read count '{clade_reads}'"
                            )
                        lineage = {r: '' for r in self.tax_ranks}
                        for _depth, rank_code, taxon_name in stack:
                            rank_name = self.KRAKEN_RANK_MAP.get(rank_code)
                            if rank_name:
                                lineage[rank_name] = taxon_name
                        counts[taxid] = counts.get(taxid, 0) + clade_count
                        tax_rows.setdefault(taxid, {'taxon_id': taxid, **lineage})

            per_sample_counts[sample_name] = counts

        if not any(per_sample_counts.values()):
            raise ValueError(
                "No species-level rows found in any kraken/bracken report; "
                "phyloseq tables cannot be built. Check the kraken2 database "
                "and bracken settings."
            )

        # Outer-join samples on taxid; taxa absent from a sample get 0, and the
        # taxonomy table covers the union of taxa across all samples.
        otu_table = pd.concat(
            [
                pd.Series(counts, name=sample, dtype=float)
                for sample, counts in per_sample_counts.items()
            ],
            axis=1
        ).fillna(0)

        tax_table = pd.DataFrame(list(tax_rows.values())).set_index('taxon_id')
        tax_table = tax_table.reindex(otu_table.index).fillna('')

        return otu_table, tax_table
    
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
            except Exception as e:
                raise ValueError(f"Failed to parse sourmash gather file {gather_file}: {e}")

            # Extract sample name from filename (<id>_smgather_<db>.with-lineages.csv)
            sample_name = gather_file.name.split('_smgather_')[0]

            # Header-only CSV: the SOURMASH module's no-match fallback
            if df.empty:
                print(
                    f"Warning: no sourmash matches for sample {sample_name} - "
                    "excluded from phyloseq tables",
                    file=sys.stderr
                )
                continue

            missing = {'name', 'unique_intersect_bp', 'scaled', 'average_abund', 'lineage'} - set(df.columns)
            if missing:
                raise ValueError(
                    f"{gather_file} is missing expected gather columns: {sorted(missing)}"
                )

            df['sample_id'] = sample_name
            all_data.append(df)

        if not all_data:
            raise ValueError(
                "No sourmash matches in any sample; phyloseq tables cannot be built. "
                "Check the sourmash database and k-mer settings (the classification "
                "report from TAXONOMY_REPORT shows per-sample classified fractions), "
                "or rerun with --taxonomic_profiler none."
            )

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
    
    # sylph-tax lineage prefixes mapped to phyloseq ranks. GTDB lineages start
    # at d__ (domain); viral taxonomies use r__ (realm) and k__ (kingdom), so
    # k__ fills Kingdom when there is no d__. t__ (the genome) is the OTU id.
    SYLPH_RANK_MAP = {
        'd': 'Kingdom', 'k': 'Kingdom', 'p': 'Phylum', 'c': 'Class',
        'o': 'Order', 'f': 'Family', 'g': 'Genus', 's': 'Species',
    }

    def parse_sylph(
        self,
        input_files: List[Path]
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Parse sylph-tax profiles plus the SYLPH_PROFILE per-sample reports.

        <id>_sylph_<db>.sylphmpa (<db> = --db-name, stripped as an exact
        suffix: a custom database label may itself contain 'sylph'):
        '#SampleID ...' line, header line, then one
        row per clade at every rank; only the genome rows (last lineage
        component t__<accession>) are used, so no abundance is counted twice.
        <id>_sylph_<db>_report.tsv: Id, Sylph DB, Unclassified, Classified,
        Read pairs.

        Args:
            input_files: .sylphmpa files and the matching *_report.tsv files

        Returns:
            Tuple of (otu_table of relative_abundance, taxonomy_table,
            read_counts): read_counts holds the estimated read pairs per
            genome, sequence_abundance / 100 x the sample's read pairs (sylph
            runs with -u, so sequence_abundance is the genome's share of all
            input reads, unknown included)
        """
        profile_suffix = f"_sylph_{self.db_name}.sylphmpa"
        profiles = sorted(f for f in input_files if f.name.endswith('.sylphmpa'))
        reports = [f for f in input_files if f.name.endswith('_report.tsv')]
        if not profiles:
            raise ValueError("No .sylphmpa files among the sylph inputs")

        read_pairs = {}
        for report_file in reports:
            df = pd.read_csv(report_file, sep='\t', dtype={'Id': str, 'Read pairs': str})
            if not {'Id', 'Read pairs'} <= set(df.columns):
                raise ValueError(f"{report_file} is missing the Id / Read pairs columns of a SYLPH_PROFILE report")
            for sample_id, pairs in zip(df['Id'], df['Read pairs']):
                # a plain non-negative integer only: '1.23457e+07' would be a
                # rounded count
                if not str(pairs).isdigit():
                    raise ValueError(f"{report_file}: Read pairs '{pairs}' for {sample_id} is not a whole number")
                read_pairs[sample_id] = int(pairs)

        required = {'clade_name', 'relative_abundance', 'sequence_abundance'}
        all_data = []
        tax_rows = {}

        for profile in profiles:
            if not profile.name.endswith(profile_suffix) or profile.name == profile_suffix:
                raise ValueError(f"{profile} is not named <sample>{profile_suffix}")
            sample_name = profile.name[:-len(profile_suffix)]
            with open(profile) as fh:
                first = fh.readline()
            if not first.startswith('#SampleID'):
                raise ValueError(f"{profile} is not a sylph-tax profile (no #SampleID line)")
            try:
                df = pd.read_csv(profile, sep='\t', skiprows=1)
            except Exception as e:
                raise ValueError(f"Failed to parse sylph-tax profile {profile}: {e}")

            missing = required - set(df.columns)
            if missing:
                raise ValueError(f"{profile} is missing expected sylph-tax columns: {sorted(missing)}")

            # Header-only profile: the SYLPH_TAXPROF no-hit fallback
            if df.empty:
                print(
                    f"Warning: no sylph genomes for sample {sample_name} - "
                    "excluded from phyloseq tables",
                    file=sys.stderr
                )
                continue

            leaves = df[df['clade_name'].str.split('|').str[-1].str.startswith('t__')].copy()
            if leaves.empty:
                raise ValueError(f"{profile} has clade rows but no t__ genome rows")
            if sample_name not in read_pairs:
                raise ValueError(f"No SYLPH_PROFILE report (Read pairs) for sample {sample_name}")

            leaves['taxon_id'] = leaves['clade_name'].str.split('|').str[-1].str[3:]
            leaves['sample_id'] = sample_name
            leaves['read_pairs'] = (
                leaves['sequence_abundance'] / 100.0 * read_pairs[sample_name]
            ).round().astype(int)
            all_data.append(leaves)

            for taxon_id, clade in leaves[['taxon_id', 'clade_name']].values:
                if taxon_id in tax_rows:
                    continue
                row = {'taxon_id': taxon_id, **{rank: '' for rank in self.tax_ranks}}
                for part in clade.split('|'):
                    prefix, sep, value = part.partition('__')
                    rank = self.SYLPH_RANK_MAP.get(prefix)
                    if sep and rank and not row[rank]:
                        row[rank] = value
                tax_rows[taxon_id] = row

        if not all_data:
            raise ValueError(
                "No sylph genomes in any sample; phyloseq tables cannot be built. "
                "Check the sylph database (the classification report from "
                "TAXONOMY_REPORT shows per-sample classified fractions), "
                "or rerun with --taxonomic_profiler none."
            )

        combined = pd.concat(all_data, ignore_index=True)
        # sylph-tax writes one t__ row per genome; a duplicate would be
        # averaged by one pivot and summed by the other
        duplicated = combined.duplicated(['taxon_id', 'sample_id'])
        if duplicated.any():
            dup = combined.loc[duplicated, ['sample_id', 'taxon_id']].iloc[0]
            raise ValueError(
                f"Genome {dup['taxon_id']} appears twice in the sylph-tax profile of sample {dup['sample_id']}"
            )
        otu_table = combined.pivot_table(
            index='taxon_id', columns='sample_id',
            values='relative_abundance', fill_value=0
        )
        read_counts = combined.pivot_table(
            index='taxon_id', columns='sample_id',
            values='read_pairs', fill_value=0, aggfunc='sum'
        ).reindex(index=otu_table.index, columns=otu_table.columns).fillna(0).astype(int)

        tax_table = pd.DataFrame(list(tax_rows.values())).set_index('taxon_id')
        tax_table = tax_table.reindex(otu_table.index).fillna('')

        return otu_table, tax_table, read_counts

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
        save_hdf5: bool = True,
        extra_tables: Optional[Dict[str, pd.DataFrame]] = None
    ) -> Dict[str, Path]:
        """
        Save phyloseq-compatible tables.
        
        Args:
            otu_table: OTU abundance table
            tax_table: Taxonomy table
            sample_metadata: Sample metadata
            save_hdf5: Whether to save HDF5 format
            extra_tables: Further taxa x samples tables, written as
                <prefix>_<name>.tsv (sylph: read_counts); the names must not
                end in otu_table, which the module's otu_table output globs
            
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

        for name, table in (extra_tables or {}).items():
            extra_path = self.output_dir / f"{prefix}_{name}.tsv"
            table.to_csv(extra_path, sep='\t')
            output_files[name] = extra_path
            print(f"Saved {name} table: {extra_path}")
        
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
            
            top_taxa = self.select_plot_taxa(rel_abundance)
            
            # Samples top to bottom in natural order (V2 before V10)
            samples = sorted(rel_abundance.columns, key=natural_key, reverse=True)
            plot_data = rel_abundance.loc[top_taxa, samples]
            
            # Add "Other" category
            other = rel_abundance.loc[~rel_abundance.index.isin(top_taxa), samples].sum(axis=0)
            plot_data.loc['Other'] = other
            
            # Create plot; height and label size follow the number of samples
            n_samples = len(samples)
            fig, ax = plt.subplots(figsize=(16, max(6, 0.3 * n_samples + 2)))
            
            # Generate color palette; "Other" is grey
            palette = "tab10" if len(top_taxa) <= 10 else "tab20"
            colors = sns.color_palette(palette, len(top_taxa)) + [(0.8, 0.8, 0.8)]
            
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
            ax.tick_params(axis='x', labelsize=18)
            ax.tick_params(axis='y', labelsize=max(6, min(18, 500 / max(n_samples, 1))))
            
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
            plot_path = self.output_dir / f"{self.profiler}_{self.db_name}_{tax_level.lower()}_bar.png"
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Created {tax_level} abundance plot: {plot_path}")
            return plot_path
            
        except Exception as e:
            print(f"Warning: Failed to create {tax_level} plot: {e}", file=sys.stderr)
            return None
    
    def select_plot_taxa(self, rel_abundance: pd.DataFrame) -> List[str]:
        """
        Choose at most top_n taxa so that every sample's dominant taxa are shown.

        Ranking by mean abundance alone hides taxa that dominate a few samples
        (they end up in 'Other'). Take the union of each sample's k most
        abundant taxa for the largest k that fits in top_n, then fill the
        remaining slots by mean abundance.
        """
        mean_abundance = rel_abundance.mean(axis=1).sort_values(ascending=False)
        chosen: List[str] = []
        for k in range(1, self.top_n + 1):
            union = set()
            for sample in rel_abundance.columns:
                col = rel_abundance[sample]
                union.update(col[col > 0].nlargest(k).index)
            if len(union) > self.top_n:
                if k == 1:
                    # more samples with distinct dominant taxa than slots:
                    # keep the dominant taxa that reach the highest abundance
                    peak = rel_abundance.loc[list(union)].max(axis=1)
                    chosen = list(peak.sort_values(ascending=False).index[:self.top_n])
                break
            chosen = [t for t in mean_abundance.index if t in union]
        for taxon in mean_abundance.index:
            if len(chosen) >= self.top_n:
                break
            if taxon not in chosen and mean_abundance[taxon] > 0:
                chosen.append(taxon)
        return chosen

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
            input_files: List of input files (kraken-style bracken reports, Sourmash
                gather, or sylph-tax profiles plus SYLPH_PROFILE reports)
            output_format: 'tables', 'hdf5', or 'both'
            
        Returns:
            Dictionary with output information
        """
        print(f"Processing {len(input_files)} input files with {self.profiler} profiler")
        
        # Parse input files based on profiler
        extra_tables = None
        if self.profiler == 'kraken2':
            # For Kraken2, we expect per-sample kraken-style bracken reports
            otu_table, tax_table = self.parse_kraken_reports(input_files)
        elif self.profiler == 'sylph':
            # For Sylph: relative abundance in the OTU table, estimated read
            # pairs per genome in a separate table
            otu_table, tax_table, read_counts = self.parse_sylph(input_files)
            extra_tables = {'read_counts': read_counts}
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
            save_hdf5=save_hdf5,
            extra_tables=extra_tables
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
  # Kraken2/Bracken kraken-style reports
  taxonomy_phyloseq.py --profiler kraken2 --input-files *_bracken.txt \\
      --db-name silva --output-dir . --format both \\
      --plot-levels Phylum,Family,Genus,Species

  # Sourmash gather CSV files
  taxonomy_phyloseq.py --profiler sourmash --input-files *_smgather*.csv \\
      --db-name gtdb --output-dir . --format both \\
      --plot-levels Phylum,Family,Genus

  # sylph-tax profiles plus the SYLPH_PROFILE per-sample reports
  taxonomy_phyloseq.py --profiler sylph --input-files *.sylphmpa *_report.tsv \\
      --db-name gtdb_r232 --output-dir . --format both \\
      --plot-levels Phylum,Family,Genus,Species
        """
    )
    
    parser.add_argument(
        '--profiler',
        required=True,
        choices=PROFILERS,
        help='Taxonomic profiler used'
    )
    
    parser.add_argument(
        '--input-files',
        required=True,
        nargs='+',
        type=Path,
        help='Input files (kraken-style bracken reports for kraken2, gather CSV for sourmash)'
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
        help='Maximum taxa shown in plots; every sample\'s most abundant taxa are included first (default: 10)'
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
