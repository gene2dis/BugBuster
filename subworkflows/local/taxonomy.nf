/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TAXONOMY PROFILING SUBWORKFLOW - OPTIMIZED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Taxonomic classification at read level using Kraken2 or Sourmash
    
    Improvements:
    - Unified Python-based reporting and phyloseq generation
    - Phyloseq-compatible TSV tables for interoperability
    - Optional R phyloseq object generation
    - Removed R container dependencies
    - Better parallelization and resource utilization
    - Reduced from 6 modules to 2 core modules (+1 optional)
----------------------------------------------------------------------------------------
*/

include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { SOURMASH             } from '../../modules/local/sourmash/main'
include { TAXONOMY_REPORT      } from '../../modules/local/taxonomy_report/main'
include { TAXONOMY_PHYLOSEQ    } from '../../modules/local/taxonomy_phyloseq/main'
include { PHYLOSEQ_CONVERTER   } from '../../modules/local/phyloseq_converter/main'

workflow TAXONOMY {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    reads_report    // channel: path(report)
    kraken_db       // channel: path(kraken_db)
    sourmash_db     // channel: path(sourmash_db)

    main:
    ch_versions = Channel.empty()

    //
    // Kraken2 taxonomic profiling
    //
    if ( params.taxonomic_profiler == "kraken2" ) {
        // Run Kraken2 classification
        KRAKEN2(
            reads,
            kraken_db.first(),
            false,  // save_output_fastqs
            false   // save_reads_assignment
        )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        
        // Run Bracken abundance estimation
        BRACKEN(
            KRAKEN2.out.report,
            kraken_db.first()
        )
        ch_versions = ch_versions.mix(BRACKEN.out.versions.first())
        
        // Generate unified taxonomy report
        TAXONOMY_REPORT(
            KRAKEN2.out.report.map { _meta, report -> report }.collect(),
            reads_report.flatten().filter { file -> file.name.endsWith('.csv') }.first(),
            'kraken2',
            params.kraken_db_used
        )
        ch_versions = ch_versions.mix(TAXONOMY_REPORT.out.versions)
        
        // Generate phyloseq-compatible tables and plots
        TAXONOMY_PHYLOSEQ(
            BRACKEN.out.txt.map { _meta, report -> report }.collect(),
            'kraken2',
            params.kraken_db_used,
            params.taxonomy_plot_levels ?: 'Phylum,Family,Genus,Species',
            params.taxonomy_top_n_taxa ?: 10
        )
        ch_versions = ch_versions.mix(TAXONOMY_PHYLOSEQ.out.versions)
        
        // Optional: Convert to R phyloseq object
        if (params.create_phyloseq_rds) {
            PHYLOSEQ_CONVERTER(
                TAXONOMY_PHYLOSEQ.out.otu_table,
                TAXONOMY_PHYLOSEQ.out.tax_table,
                TAXONOMY_PHYLOSEQ.out.sample_metadata,
                params.kraken_db_used
            )
            ch_versions = ch_versions.mix(PHYLOSEQ_CONVERTER.out.versions)
        }
    }

    //
    // Sourmash taxonomic profiling
    //
    if ( params.taxonomic_profiler == "sourmash" ) {
        // Run Sourmash classification
        ch_sm_taxonomy = SOURMASH(
            reads.combine(sourmash_db),
            params.sourmash_db_name,
            params.sourmash_tax_rank
        )
        
        // Generate unified taxonomy report
        TAXONOMY_REPORT(
            ch_sm_taxonomy.report.collect(),
            reads_report.flatten().filter { file -> file.name.endsWith('.csv') }.first(),
            'sourmash',
            params.sourmash_db_name
        )
        ch_versions = ch_versions.mix(TAXONOMY_REPORT.out.versions)
        
        // Generate phyloseq-compatible tables and plots
        TAXONOMY_PHYLOSEQ(
            ch_sm_taxonomy.sourmash_gather.collect(),
            'sourmash',
            params.sourmash_db_name,
            params.taxonomy_plot_levels ?: 'Phylum,Family,Genus,Species',
            params.taxonomy_top_n_taxa ?: 10
        )
        ch_versions = ch_versions.mix(TAXONOMY_PHYLOSEQ.out.versions)
        
        // Optional: Convert to R phyloseq object
        if (params.create_phyloseq_rds) {
            PHYLOSEQ_CONVERTER(
                TAXONOMY_PHYLOSEQ.out.otu_table,
                TAXONOMY_PHYLOSEQ.out.tax_table,
                TAXONOMY_PHYLOSEQ.out.sample_metadata,
                params.sourmash_db_name
            )
            ch_versions = ch_versions.mix(PHYLOSEQ_CONVERTER.out.versions)
        }
    }

    emit:
    versions = ch_versions // channel: path(versions.yml)
}
