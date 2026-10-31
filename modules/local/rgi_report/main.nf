process RGI_REPORT {
    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"
    
    tag "rgi_report"
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
        // --rgi_summary_min_coverage: a gene counts in the *_cov<N> columns and the
        // plots when reads cover at least N% of its best-covered reference sequence
        // (rgi bwt itself lists every gene with a mapped read); the raw RGI tables
        // and the unfiltered columns are unchanged
        // != null, not ?: (Groovy reads 0 as false, and 0 is a valid cutoff)
        def min_cov = params.rgi_summary_min_coverage != null ? params.rgi_summary_min_coverage : 80
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
        # RGI_KMER's k-mer JSON (rgi kmer_query --bwt writes no *mer_analysis.txt);
        # counted for the log only, the summary and plots use the gene tables
        kmer_files = [f for f in Path('.').glob('*mer_analysis.json')]
        
        print(f"Found {len(allele_files)} allele files")
        print(f"Found {len(gene_files)} gene files")
        print(f"Found {len(kmer_files)} k-mer files")

        def aro_key(value):
            # ARO accessions as text; a table with NA rows reads them as floats (3000191.0)
            try:
                return str(int(float(value)))
            except (TypeError, ValueError):
                return ''

        # Best coverage of each ARO per sample, from the allele tables
        min_cov = float("${min_cov}")
        cov_tag = f"cov{min_cov:g}"
        detected = {}
        for allele_file in allele_files:
            sample_id = allele_file.name.replace('_rgi_bwt.allele_mapping_data.txt', '')
            try:
                a = pd.read_csv(allele_file, sep='\t')
            except Exception as e:
                print(f"Warning: Could not process {allele_file}: {e}")
                continue
            if {'ARO Accession', 'Percent Coverage'} <= set(a.columns):
                cov = pd.to_numeric(a['Percent Coverage'], errors='coerce')
                best = cov.groupby(a['ARO Accession'].map(aro_key)).max()
                detected[sample_id] = set(best.index[best >= min_cov])
        
        # Process gene-level data. A sample where rgi bwt mapped too few reads
        # gets a placeholder table with one all-NA row (rgi_bwt_lowcov_fallback.sh);
        # that row is dropped, so the sample is reported with zero genes
        samples = []
        gene_data_list = []
        for gene_file in gene_files:
            try:
                sample_id = gene_file.name.replace('_rgi_bwt.gene_mapping_data.txt', '')
                df = pd.read_csv(gene_file, sep='\\t')
                if 'ARO Term' in df.columns:
                    df = df[df['ARO Term'].notna()].copy()
                df['Sample'] = sample_id
                samples.append(sample_id)
                gene_data_list.append(df)
            except Exception as e:
                print(f"Warning: Could not process {gene_file}: {e}")
        samples = sorted(samples)
        
        summary_columns = ['Sample', 'Total_AMR_Genes', 'Total_Mapped_Reads',
                           'Unique_Drug_Classes', 'Unique_AMR_Families',
                           f'AMR_Genes_{cov_tag}', f'Mapped_Reads_{cov_tag}',
                           f'Unique_Drug_Classes_{cov_tag}', f'Unique_AMR_Families_{cov_tag}']
        if gene_data_list:
            combined_gene_data = pd.concat(gene_data_list, ignore_index=True)
            if 'ARO Accession' in combined_gene_data.columns:
                keep = [aro_key(aro) in detected.get(s, set()) for s, aro in
                        zip(combined_gene_data['Sample'], combined_gene_data['ARO Accession'])]
            else:
                keep = [False] * len(combined_gene_data)
            filtered = combined_gene_data[pd.Series(keep, index=combined_gene_data.index, dtype=bool)]
            
            # Create summary report, one row per sample (zero-hit samples included)
            summary_data = []
            for sample in samples:
                sample_data = combined_gene_data[combined_gene_data['Sample'] == sample]
                summary_data.append({
                    'Sample': sample,
                    'Total_AMR_Genes': len(sample_data),
                    'Total_Mapped_Reads': int(sample_data['All Mapped Reads'].sum()) if 'All Mapped Reads' in sample_data.columns else 0,
                    'Unique_Drug_Classes': sample_data['Drug Class'].nunique() if 'Drug Class' in sample_data.columns else 0,
                    'Unique_AMR_Families': sample_data['AMR Gene Family'].nunique() if 'AMR Gene Family' in sample_data.columns else 0
                })
                fsd = filtered[filtered['Sample'] == sample]
                summary_data[-1].update({
                    f'AMR_Genes_{cov_tag}': len(fsd),
                    f'Mapped_Reads_{cov_tag}': int(fsd['All Mapped Reads'].sum()) if 'All Mapped Reads' in fsd.columns else 0,
                    f'Unique_Drug_Classes_{cov_tag}': fsd['Drug Class'].nunique() if 'Drug Class' in fsd.columns else 0,
                    f'Unique_AMR_Families_{cov_tag}': fsd['AMR Gene Family'].nunique() if 'AMR Gene Family' in fsd.columns else 0,
                })
            
            summary_df = pd.DataFrame(summary_data, columns=summary_columns)
            summary_df.to_csv('RGI_summary_report.csv', index=False)
            
            # Generate plots from the genes at or above the coverage cutoff
            combined_gene_data = filtered
            cut = f" (>= {min_cov:g}% gene coverage)"
            try:
                # Plot 1: AMR Gene Family distribution
                if 'AMR Gene Family' in combined_gene_data.columns and combined_gene_data['AMR Gene Family'].notna().any():
                    top_families = combined_gene_data['AMR Gene Family'].value_counts().head(15)
                    plt.figure(figsize=(14, 8))
                    top_families.plot(kind='barh')
                    plt.xlabel('Number of Occurrences')
                    plt.ylabel('AMR Gene Family')
                    plt.title('Top 15 AMR Gene Families Detected' + cut)
                    plt.tight_layout()
                    plt.savefig('RGI_amr_gene_family_distribution.png', dpi=300)
                    plt.close()
                
                # Plot 2: Drug Class profile
                if 'Drug Class' in combined_gene_data.columns and combined_gene_data['Drug Class'].notna().any():
                    drug_classes = combined_gene_data['Drug Class'].str.split(';').explode()
                    top_drugs = drug_classes.value_counts().head(15)
                    plt.figure(figsize=(14, 8))
                    top_drugs.plot(kind='barh', color='coral')
                    plt.xlabel('Number of Occurrences')
                    plt.ylabel('Drug Class')
                    plt.title('Top 15 Drug Classes with Resistance' + cut)
                    plt.tight_layout()
                    plt.savefig('RGI_drug_class_profile.png', dpi=300)
                    plt.close()
                
                # Plot 3: Resistance Mechanism
                if 'Resistance Mechanism' in combined_gene_data.columns and combined_gene_data['Resistance Mechanism'].notna().any():
                    mechanisms = combined_gene_data['Resistance Mechanism'].str.split(';').explode()
                    top_mechanisms = mechanisms.value_counts().head(10)
                    plt.figure(figsize=(12, 8))
                    top_mechanisms.plot(kind='barh', color='steelblue')
                    plt.xlabel('Number of Occurrences')
                    plt.ylabel('Resistance Mechanism')
                    plt.title('Top 10 Resistance Mechanisms' + cut)
                    plt.tight_layout()
                    plt.savefig('RGI_resistance_mechanisms.png', dpi=300)
                    plt.close()
                
                # Plot 4: Per-sample AMR gene count (zero-hit samples shown as 0)
                plt.figure(figsize=(12, 6))
                sample_counts = combined_gene_data.groupby('Sample').size().reindex(samples, fill_value=0)
                sample_counts.plot(kind='bar', color='mediumseagreen')
                plt.xlabel('Sample')
                plt.ylabel('Number of AMR Genes')
                plt.title('AMR Genes Detected per Sample' + cut)
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig('RGI_sample_amr_counts.png', dpi=300)
                plt.close()
                
            except Exception as e:
                print(f"Warning: Could not generate plots: {e}")
        else:
            # Create empty summary if no data
            summary_df = pd.DataFrame(columns=summary_columns)
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
            python: 3.13.18
            pandas: 3.1.0
        END_VERSIONS
        """
}
