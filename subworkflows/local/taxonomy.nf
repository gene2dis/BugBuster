/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TAXONOMY PROFILING SUBWORKFLOW - OPTIMIZED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Taxonomic classification at read level using Kraken2, Sourmash or Sylph
    
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
include { SYLPH_PROFILE        } from '../../modules/local/sylph_profile/main'
include { SYLPH_TAXPROF        } from '../../modules/local/sylph_taxprof/main'
include { TAXONOMY_REPORT      } from '../../modules/local/taxonomy_report/main'
include { TAXONOMY_PHYLOSEQ    } from '../../modules/local/taxonomy_phyloseq/main'
include { PHYLOSEQ_CONVERTER   } from '../../modules/local/phyloseq_converter/main'
include { flagOn; sylphDbLabel; kraken2DbLabel; sourmashDbLabel } from './utils_params'

workflow TAXONOMY {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    reads_report    // channel: path(report)
    kraken_db       // channel: path(kraken_db)
    sourmash_db     // channel: path(sourmash_db)
    sylph_db        // channel: path(sylph_db)
    sylph_taxonomy  // channel: path(sylph-tax metadata)

    main:
    ch_versions = channel.empty()
    // For the results summary
    ch_tax_reads_report = channel.empty()
    ch_tax_tables       = channel.empty()
    ch_sylph_profiles   = channel.empty()

    //
    // Kraken2 taxonomic profiling
    //
    if ( params.taxonomic_profiler == "kraken2" ) {
        // Run Kraken2 classification on the paired reads only: the nf-core
        // module passes the whole read list after --paired, which requires an
        // even number of files (a staged singleton would abort kraken2)
        KRAKEN2(
            reads.map { meta, reads_files ->
                [ meta, reads_files.size() > 2 ? reads_files[0..1] : reads_files ]
            },
            kraken_db.first(),
            false,  // save_output_fastqs
            false   // save_reads_assignment
        )
        // KRAKEN2 and BRACKEN report their versions on the versions topic
        // (collected in main.nf), not through ch_versions

        // Run Bracken abundance estimation
        BRACKEN(
            KRAKEN2.out.report,
            kraken_db.first()
        )
        
        // Generate unified taxonomy report
        TAXONOMY_REPORT(
            KRAKEN2.out.report.map { _meta, report -> report }.collect(sort: true),
            reads_report.flatten().filter { file -> file.name.endsWith('.csv') }.first(),
            'kraken2',
            kraken2DbLabel()
        )
        ch_versions = ch_versions.mix(TAXONOMY_REPORT.out.versions)
        
        // Generate phyloseq-compatible tables and plots
        TAXONOMY_PHYLOSEQ(
            BRACKEN.out.txt.map { _meta, report -> report }.collect(sort: true),
            'kraken2',
            kraken2DbLabel(),
            params.taxonomy_plot_levels ?: 'Phylum,Family,Genus,Species',
            params.taxonomy_top_n_taxa ?: 10
        )
        ch_versions = ch_versions.mix(TAXONOMY_PHYLOSEQ.out.versions)
        
        // Optional: Convert to R phyloseq object
        if (flagOn(params.create_phyloseq_rds)) {
            PHYLOSEQ_CONVERTER(
                TAXONOMY_PHYLOSEQ.out.otu_table,
                TAXONOMY_PHYLOSEQ.out.tax_table,
                TAXONOMY_PHYLOSEQ.out.sample_metadata,
                kraken2DbLabel()
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
            sourmashDbLabel(),
            params.sourmash_tax_rank
        )
        ch_versions = ch_versions.mix(ch_sm_taxonomy.versions.first())

        // Generate unified taxonomy report
        TAXONOMY_REPORT(
            ch_sm_taxonomy.report.collect(sort: true),
            reads_report.flatten().filter { file -> file.name.endsWith('.csv') }.first(),
            'sourmash',
            sourmashDbLabel()
        )
        ch_versions = ch_versions.mix(TAXONOMY_REPORT.out.versions)
        
        // Generate phyloseq-compatible tables and plots
        TAXONOMY_PHYLOSEQ(
            ch_sm_taxonomy.sourmash_gather.collect(sort: true),
            'sourmash',
            sourmashDbLabel(),
            params.taxonomy_plot_levels ?: 'Phylum,Family,Genus,Species',
            params.taxonomy_top_n_taxa ?: 10
        )
        ch_versions = ch_versions.mix(TAXONOMY_PHYLOSEQ.out.versions)
        
        // Optional: Convert to R phyloseq object
        if (flagOn(params.create_phyloseq_rds)) {
            PHYLOSEQ_CONVERTER(
                TAXONOMY_PHYLOSEQ.out.otu_table,
                TAXONOMY_PHYLOSEQ.out.tax_table,
                TAXONOMY_PHYLOSEQ.out.sample_metadata,
                sourmashDbLabel()
            )
            ch_versions = ch_versions.mix(PHYLOSEQ_CONVERTER.out.versions)
        }
    }

    //
    // Sylph taxonomic profiling
    //
    if ( params.taxonomic_profiler == "sylph" ) {
        // Database label for file names, the Sylph DB column and plots
        def sylph_db_name = sylphDbLabel()

        // Paired reads only, as for kraken2: sylph profiles every -r file as
        // a separate sample, so a singleton cannot join its pair
        SYLPH_PROFILE(
            reads.map { meta, reads_files ->
                [ meta, reads_files.size() > 2 ? reads_files[0..1] : reads_files ]
            },
            sylph_db.first(),
            sylph_db_name
        )
        ch_versions = ch_versions.mix(SYLPH_PROFILE.out.versions.first())

        // Genome-level profile -> MetaPhlAn-style lineage profile
        SYLPH_TAXPROF(
            SYLPH_PROFILE.out.profile,
            sylph_taxonomy.first(),
            sylph_db_name
        )
        ch_versions = ch_versions.mix(SYLPH_TAXPROF.out.versions.first())

        // Generate unified taxonomy report
        TAXONOMY_REPORT(
            SYLPH_PROFILE.out.report.collect(sort: true),
            reads_report.flatten().filter { file -> file.name.endsWith('.csv') }.first(),
            'sylph',
            sylph_db_name
        )
        ch_versions = ch_versions.mix(TAXONOMY_REPORT.out.versions)

        // Phyloseq-compatible tables and plots: relative abundance in the OTU
        // table, plus estimated read pairs per genome, which need each
        // sample's read-pair count from its SYLPH_PROFILE report
        TAXONOMY_PHYLOSEQ(
            SYLPH_TAXPROF.out.sylphmpa.map { _meta, profile -> profile }
                .mix(SYLPH_PROFILE.out.report)
                .collect(sort: true),
            'sylph',
            sylph_db_name,
            params.taxonomy_plot_levels ?: 'Phylum,Family,Genus,Species',
            params.taxonomy_top_n_taxa ?: 10
        )
        ch_versions = ch_versions.mix(TAXONOMY_PHYLOSEQ.out.versions)

        // Optional: Convert to R phyloseq object
        if (flagOn(params.create_phyloseq_rds)) {
            PHYLOSEQ_CONVERTER(
                TAXONOMY_PHYLOSEQ.out.otu_table,
                TAXONOMY_PHYLOSEQ.out.tax_table,
                TAXONOMY_PHYLOSEQ.out.sample_metadata,
                sylph_db_name
            )
            ch_versions = ch_versions.mix(PHYLOSEQ_CONVERTER.out.versions)
        }
    }

    // TAXONOMY_REPORT runs in exactly one branch above, so its outputs are
    // read once here rather than in every branch
    if ( params.taxonomic_profiler in ['kraken2', 'sourmash', 'sylph'] ) {
        ch_tax_reads_report = TAXONOMY_REPORT.out.report
        ch_tax_tables = TAXONOMY_PHYLOSEQ.out.otu_table
            .mix(TAXONOMY_PHYLOSEQ.out.tax_table, TAXONOMY_PHYLOSEQ.out.read_counts)
    }
    if ( params.taxonomic_profiler == "sylph" ) {
        ch_sylph_profiles = SYLPH_PROFILE.out.profile.map { _meta, profile -> profile }
    }

    emit:
    reads_report   = ch_tax_reads_report // channel: path(Reads_report.csv) with the classified fraction
    tables         = ch_tax_tables       // channel: path(<label>_{otu_table,tax_table,read_counts}.tsv)
    sylph_profiles = ch_sylph_profiles   // channel: path(<sample>_sylph_<db>.tsv)
    versions       = ch_versions         // channel: path(versions.yml)
}
