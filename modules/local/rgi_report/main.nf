process RGI_REPORT {
    container 'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0'
    
    label 'process_low'

    input:
        path(allele_files)
        path(gene_files)
        path(kmer_files)

    output:
        path("RGI_summary_report.csv"), emit: summary
        path("RGI_*.png"), emit: plots, optional: true
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        #!/usr/bin/env python3
        
        import pandas as pd
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        from pathlib import Path
        import json
        import sys
        
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.figsize'] = (12, 8)
        
        # Collect all allele mapping files
        allele_files = [f for f in Path('.').glob('*.allele_mapping_data.txt')]
        gene_files = [f for f in Path('.').glob('*.gene_mapping_data.txt')]
        kmer_files = [f for f in Path('.').glob('*_61mer_analysis.txt')]
        
        print(f"Found {len(allele_files)} allele files")
        print(f"Found {len(gene_files)} gene files")
        print(f"Found {len(kmer_files)} k-mer files")
        
        # Process gene-level data
        gene_data_list = []
        for gene_file in gene_files:
            try:
                sample_id = gene_file.name.replace('_rgi_bwt.gene_mapping_data.txt', '')
                df = pd.read_csv(gene_file, sep='\\t')
                df['Sample'] = sample_id
                gene_data_list.append(df)
            except Exception as e:
                print(f"Warning: Could not process {gene_file}: {e}")
        
        if gene_data_list:
            combined_gene_data = pd.concat(gene_data_list, ignore_index=True)
            
            # Create summary report
            summary_data = []
            for sample in combined_gene_data['Sample'].unique():
                sample_data = combined_gene_data[combined_gene_data['Sample'] == sample]
                summary_data.append({
                    'Sample': sample,
                    'Total_AMR_Genes': len(sample_data),
                    'Total_Mapped_Reads': sample_data['All Mapped Reads'].sum() if 'All Mapped Reads' in sample_data.columns else 0,
                    'Unique_Drug_Classes': sample_data['Drug Class'].nunique() if 'Drug Class' in sample_data.columns else 0,
                    'Unique_AMR_Families': sample_data['AMR Gene Family'].nunique() if 'AMR Gene Family' in sample_data.columns else 0
                })
            
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv('RGI_summary_report.csv', index=False)
            
            # Generate plots
            try:
                # Plot 1: AMR Gene Family distribution
                if 'AMR Gene Family' in combined_gene_data.columns:
                    top_families = combined_gene_data['AMR Gene Family'].value_counts().head(15)
                    plt.figure(figsize=(14, 8))
                    top_families.plot(kind='barh')
                    plt.xlabel('Number of Occurrences')
                    plt.ylabel('AMR Gene Family')
                    plt.title('Top 15 AMR Gene Families Detected')
                    plt.tight_layout()
                    plt.savefig('RGI_amr_gene_family_distribution.png', dpi=300)
                    plt.close()
                
                # Plot 2: Drug Class profile
                if 'Drug Class' in combined_gene_data.columns:
                    drug_classes = combined_gene_data['Drug Class'].str.split(';').explode()
                    top_drugs = drug_classes.value_counts().head(15)
                    plt.figure(figsize=(14, 8))
                    top_drugs.plot(kind='barh', color='coral')
                    plt.xlabel('Number of Occurrences')
                    plt.ylabel('Drug Class')
                    plt.title('Top 15 Drug Classes with Resistance')
                    plt.tight_layout()
                    plt.savefig('RGI_drug_class_profile.png', dpi=300)
                    plt.close()
                
                # Plot 3: Resistance Mechanism
                if 'Resistance Mechanism' in combined_gene_data.columns:
                    mechanisms = combined_gene_data['Resistance Mechanism'].str.split(';').explode()
                    top_mechanisms = mechanisms.value_counts().head(10)
                    plt.figure(figsize=(12, 8))
                    top_mechanisms.plot(kind='barh', color='steelblue')
                    plt.xlabel('Number of Occurrences')
                    plt.ylabel('Resistance Mechanism')
                    plt.title('Top 10 Resistance Mechanisms')
                    plt.tight_layout()
                    plt.savefig('RGI_resistance_mechanisms.png', dpi=300)
                    plt.close()
                
                # Plot 4: Per-sample AMR gene count
                plt.figure(figsize=(12, 6))
                sample_counts = combined_gene_data.groupby('Sample').size()
                sample_counts.plot(kind='bar', color='mediumseagreen')
                plt.xlabel('Sample')
                plt.ylabel('Number of AMR Genes')
                plt.title('AMR Genes Detected per Sample')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig('RGI_sample_amr_counts.png', dpi=300)
                plt.close()
                
            except Exception as e:
                print(f"Warning: Could not generate plots: {e}")
        else:
            # Create empty summary if no data
            summary_df = pd.DataFrame(columns=['Sample', 'Total_AMR_Genes', 'Total_Mapped_Reads', 
                                              'Unique_Drug_Classes', 'Unique_AMR_Families'])
            summary_df.to_csv('RGI_summary_report.csv', index=False)
        
        # Create versions file
        with open('versions.yml', 'w') as f:
            f.write('"${task.process}":\\n')
            f.write(f'    python: {sys.version.split()[0]}\\n')
            f.write(f'    pandas: {pd.__version__}\\n')
            f.write(f'    matplotlib: {matplotlib.__version__}\\n')
        
        print("RGI report generation completed")
        """

    stub:
        """
        touch RGI_summary_report.csv
        touch RGI_amr_gene_family_distribution.png
        touch RGI_drug_class_profile.png
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.9.0
            pandas: 1.3.0
        END_VERSIONS
        """
}
