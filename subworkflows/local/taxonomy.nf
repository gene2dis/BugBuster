/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TAXONOMY PROFILING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Taxonomic classification at read level using Kraken2 or Sourmash
----------------------------------------------------------------------------------------
*/

include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { KRAKEN_BIOM          } from '../../modules/local/kraken_biom/main'
include { KRAKEN_TO_PHYLOSEQ   } from '../../modules/local/kraken_to_phyloseq/main'
include { TAX_REPORT_KRAKEN2   } from '../../modules/local/tax_report_kraken2/main'
include { SOURMASH             } from '../../modules/local/sourmash/main'
include { SOURMASH_TO_PHYLOSEQ } from '../../modules/local/sourmash_to_phyloseq/main'
include { TAX_REPORT_SOURMASH  } from '../../modules/local/tax_report_sourmash/main'

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
        //
        // Run nf-core Kraken2 classification
        // nf-core KRAKEN2_KRAKEN2 signature:
        //   input:  tuple val(meta), path(reads) + path(db) + val save_output_fastqs + val save_reads_assignment
        //   output: tuple val(meta), path('*report.txt'), emit: report
        //
        KRAKEN2(
            reads,
            kraken_db.first(),
            false,  // save_output_fastqs
            false   // save_reads_assignment
        )
        
        // Collect versions
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
        
        //
        // Generate Kraken2 summary report (replaces local module's TSV generation)
        // Extract classification stats from kraken report for TAX_REPORT_KRAKEN2
        //
        ch_kraken_report = KRAKEN2.out.report
        
        // Generate taxonomy report
        TAX_REPORT_KRAKEN2(
            ch_kraken_report
                .map { _meta, report -> report }
                .concat(reads_report)
                .collect()
        )

        //
        // Bracken abundance estimation with nf-core module
        // nf-core BRACKEN_BRACKEN signature:
        //   input:  tuple val(meta), path(kraken_report) + path database
        //   output: tuple val(meta), path(bracken_report), emit: reports
        //
        BRACKEN(
            ch_kraken_report,
            kraken_db.first()
        )
        
        // Collect versions
        ch_versions = ch_versions.mix(BRACKEN.out.versions.first())
        
        // Extract bracken reports for KRAKEN_BIOM (expects path only)
        ch_bracken_reports = BRACKEN.out.txt
            .map { _meta, report -> report }

        // Generate BIOM file
        ch_kraken_biom = KRAKEN_BIOM(
            ch_bracken_reports.collect(),
            params.kraken_db_used
        )

        // Generate Phyloseq object and abundance plots
        KRAKEN_TO_PHYLOSEQ(ch_kraken_biom)
    }

    //
    // Sourmash taxonomic profiling
    //
    if ( params.taxonomic_profiler == "sourmash" ) {
        // Run Sourmash classification
        ch_sm_taxonomy = SOURMASH(
            reads.combine(sourmash_db.collect()),
            params.sourmash_db_name,
            params.sourmash_tax_rank
        )

        // Generate taxonomy report
        TAX_REPORT_SOURMASH(
            ch_sm_taxonomy.report
                .concat(reads_report)
                .collect()
        )

        // Generate Phyloseq object
        SOURMASH_TO_PHYLOSEQ(
            ch_sm_taxonomy.sourmash_gather.collect(),
            params.sourmash_db_name
        )
    }

    emit:
    versions = ch_versions // channel: path(versions.yml)
}
