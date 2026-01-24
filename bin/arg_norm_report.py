#!/usr/bin/env python3
"""
ARG Normalization Report Generator
Processes KARGA, KARGVA, and ARGs-OAP results to normalize ARG abundance.
Python migration of Read_arg_norm.R
"""

import pandas as pd
import glob
import sys
import os


def parse_karga_gene_idx(gene_idx):
    """Parse KARGA GeneIdx pipe-delimited string."""
    if pd.isna(gene_idx):
        return {'Compound': None, 'Meg_drug': None, 'Other_Drug': None, 'Gen': None}
    
    parts = gene_idx.split('|')
    result = {
        'Compound': parts[1] if len(parts) > 1 else None,
        'Meg_drug': parts[2] if len(parts) > 2 else None,
        'Other_Drug': parts[3] if len(parts) > 3 else None,
        'Gen': parts[4] if len(parts) > 4 else None
    }
    return result


def parse_kargva_gene_idx(gene_idx):
    """Parse KARGVA GeneIdx pipe-delimited string."""
    if pd.isna(gene_idx):
        return {
            'Mutation': None, 'Meg_drug': None, 'Drug': None, 'Other_Drug': None,
            'Other_Drug_2': None, 'Organism_1': None, 'Organism_2': None,
            'Gen_1': None, 'Gen_2': None, 'Gen_3': None
        }
    
    parts = gene_idx.split('|')
    result = {
        'Mutation': parts[1] if len(parts) > 1 else None,
        'Meg_drug': parts[6] if len(parts) > 6 else None,
        'Drug': parts[5] if len(parts) > 5 else None,
        'Other_Drug': parts[14] if len(parts) > 14 else None,
        'Other_Drug_2': parts[7] if len(parts) > 7 else None,
        'Organism_1': parts[3] if len(parts) > 3 else None,
        'Organism_2': parts[16] if len(parts) > 16 else None,
        'Gen_1': parts[8] if len(parts) > 8 else None,
        'Gen_2': parts[13].split('_')[0] if len(parts) > 13 else None,
        'Gen_3': parts[3].split(' ')[2] if len(parts) > 3 and len(parts[3].split(' ')) > 2 else None
    }
    return result


def process_karga_data(karga_df):
    """Process and filter KARGA data."""
    if karga_df.empty or karga_df['GeneIdx'].isna().all():
        return karga_df
    
    # Parse GeneIdx
    parsed = karga_df['GeneIdx'].apply(parse_karga_gene_idx)
    parsed_df = pd.DataFrame(parsed.tolist())
    karga_df = pd.concat([karga_df, parsed_df], axis=1)
    
    # Convert numeric columns
    karga_df['AverageKMerDepth'] = pd.to_numeric(karga_df['AverageKMerDepth'], errors='coerce')
    karga_df['PercentGeneCovered'] = karga_df['PercentGeneCovered'].str.rstrip('%').astype(float)
    
    # Filter by coverage
    karga_df = karga_df[karga_df['PercentGeneCovered'] >= 90].copy()
    
    # Return early if no data passes filter
    if karga_df.empty:
        return karga_df
    
    # Clean up drug names
    karga_df['Meg_drug'] = karga_df.apply(
        lambda row: 'Metal resistance' if row['Compound'] == 'Metals' else row['Meg_drug'], axis=1
    )
    karga_df['Meg_drug'] = karga_df.apply(
        lambda row: 'Biocide resistance' if row['Compound'] == 'Biocides' else row['Meg_drug'], axis=1
    )
    karga_df['Meg_drug'] = karga_df.apply(
        lambda row: 'Multi compound resistance' if row['Compound'] == 'Multi-compound' else row['Meg_drug'], axis=1
    )
    
    # Capitalize first letter
    karga_df['Meg_drug'] = karga_df['Meg_drug'].str.strip().str.capitalize()
    karga_df['Meg_drug'] = karga_df['Meg_drug'].str.replace('_', ' ')
    
    # Select final columns
    karga_df = karga_df[['sample', 'GeneIdx', 'Meg_drug', 'Compound', 'Gen', 
                         'PercentGeneCovered', 'AverageKMerDepth']]
    
    return karga_df


def process_kargva_data(kargva_df):
    """Process and filter KARGVA data."""
    if kargva_df.empty or kargva_df['GeneIdx'].isna().all():
        return kargva_df
    
    # Parse GeneIdx
    parsed = kargva_df['GeneIdx'].apply(parse_kargva_gene_idx)
    parsed_df = pd.DataFrame(parsed.tolist())
    kargva_df = pd.concat([kargva_df, parsed_df], axis=1)
    
    # Convert numeric columns
    kargva_df['KmerSNPHits'] = kargva_df['KmerSNPHits'].str.split('/').str[0].astype(float)
    kargva_df['PercentGeneCovered'] = kargva_df['PercentGeneCovered'].str.rstrip('%').astype(float)
    
    # Filter
    kargva_df = kargva_df[kargva_df['KmerSNPHits'] >= 2].copy()
    kargva_df = kargva_df[kargva_df['PercentGeneCovered'] >= 80].copy()
    
    # Return early if no data passes filter
    if kargva_df.empty:
        return kargva_df
    
    # Process organism fields
    kargva_df['Organism_1'] = kargva_df['Organism_1'].str.split('resistance to').str[1]
    kargva_df['Organism_2'] = kargva_df['Organism_2'].str.split('_resistant').str[0]
    
    # Determine Gen_fusion
    kargva_df['Gen_fusion'] = kargva_df['Gen_3'].fillna(
        kargva_df['Gen_2'].fillna(kargva_df['Gen_1'])
    )
    
    # Determine All_drug
    kargva_df['All_drug'] = kargva_df['Organism_1'].fillna(
        kargva_df.apply(lambda row: row['Drug'] if row['Other_Drug'] == 'NA' else row['Other_Drug'], axis=1)
    )
    
    # Clean up drug names
    drug_replacements = {
        'Multi-compound': 'MULTIDRUG',
        'Drugs': 'MULTIDRUG',
        ' rifampicin': 'RIFAMYCIN',
        ' fusidic acid': 'FUSIDIC_ACID',
        ' sulfonamides': 'SULFONAMIDE'
    }
    kargva_df['All_drug'] = kargva_df['All_drug'].replace(drug_replacements)
    
    # Capitalize first letter
    kargva_df['Drug_resistance'] = kargva_df['All_drug'].str.strip().str.capitalize()
    
    # Select final columns
    kargva_df = kargva_df[['sample', 'GeneIdx', 'Drug_resistance', 'Gen_fusion', 
                           'PercentGeneCovered', 'AverageKMerDepth']]
    kargva_df = kargva_df.rename(columns={'Gen_fusion': 'Gen'})
    
    return kargva_df


def main():
    """Main function to generate ARG normalization report."""
    
    # Read KARGA files
    karga_files = glob.glob('*KARGA_mappedGenes.csv')
    if not karga_files:
        print("ERROR: No KARGA files found (*KARGA_mappedGenes.csv)", file=sys.stderr)
        sys.exit(1)
    
    karga_dfs = []
    for file in karga_files:
        sample_name = file.split('_all_reads_')[0]
        df = pd.read_csv(file)
        df['sample'] = sample_name
        karga_dfs.append(df)
    
    all_karga = pd.concat(karga_dfs, ignore_index=True)
    
    # Read KARGVA files
    kargva_files = glob.glob('*KARGVA_mappedGenes.csv')
    if not kargva_files:
        print("ERROR: No KARGVA files found (*KARGVA_mappedGenes.csv)", file=sys.stderr)
        sys.exit(1)
    
    kargva_dfs = []
    for file in kargva_files:
        sample_name = file.split('_all_reads_')[0]
        df = pd.read_csv(file)
        df['sample'] = sample_name
        kargva_dfs.append(df)
    
    all_kargva = pd.concat(kargva_dfs, ignore_index=True)
    
    # Read ARGs-OAP cell counts
    args_oap_dirs = [d for d in glob.glob('*_args_oap_s1_out') if os.path.isdir(d)]
    if not args_oap_dirs:
        print("ERROR: No ARGs-OAP output directories found (*_args_oap_s1_out)", file=sys.stderr)
        sys.exit(1)
    
    cell_dfs = []
    for dir_path in args_oap_dirs:
        metadata_file = os.path.join(dir_path, 'metadata.txt')
        if os.path.exists(metadata_file):
            df = pd.read_csv(metadata_file, sep='\t')
            df['sample'] = df['sample'].str.split('_tmp_').str[0]
            cell_summary = df.groupby('sample').agg({
                'nRead': 'sum',
                'n16S': 'sum',
                'nCell': 'sum'
            }).reset_index()
            cell_dfs.append(cell_summary)
    
    if not cell_dfs:
        print("ERROR: No metadata.txt files found in ARGs-OAP directories", file=sys.stderr)
        sys.exit(1)
    
    reads_and_cells = pd.concat(cell_dfs, ignore_index=True)
    
    # Process KARGA data
    karga_filtered = process_karga_data(all_karga)
    
    # Process KARGVA data
    kargva_filtered = process_kargva_data(all_kargva)
    
    # Join with read/cell counts
    karga_with_reads = reads_and_cells.merge(karga_filtered, on='sample', how='left')
    kargva_with_reads = reads_and_cells.merge(kargva_filtered, on='sample', how='left')
    
    # Normalize KARGA
    if not karga_filtered.empty and not karga_filtered['GeneIdx'].isna().all():
        karga_norm = karga_with_reads.copy()
        karga_norm['CPM_ARGs'] = karga_norm['AverageKMerDepth'] * (1 / (karga_norm['nRead'] / 1e6))
        karga_norm['copies_per_cell'] = karga_norm['AverageKMerDepth'] / karga_norm['nCell']
    else:
        karga_norm = karga_with_reads
    
    # Normalize KARGVA
    if not kargva_filtered.empty and not kargva_filtered['GeneIdx'].isna().all():
        kargva_norm = kargva_with_reads.copy()
        kargva_norm['CPM_ARGs'] = kargva_norm['AverageKMerDepth'] * (1 / (kargva_norm['nRead'] / 1e6))
        kargva_norm['copies_per_cell'] = kargva_norm['AverageKMerDepth'] / kargva_norm['nCell']
    else:
        kargva_norm = kargva_with_reads
    
    # Write outputs
    karga_norm.to_csv('KARGA_norm.csv', index=False)
    kargva_norm.to_csv('KARGVA_norm.csv', index=False)
    
    print(f"Successfully generated KARGA_norm.csv ({len(karga_norm)} rows)")
    print(f"Successfully generated KARGVA_norm.csv ({len(kargva_norm)} rows)")


if __name__ == '__main__':
    main()
