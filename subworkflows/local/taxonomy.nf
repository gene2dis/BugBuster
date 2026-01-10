/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TAXONOMY PROFILING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Taxonomic classification at read level using Kraken2 or Sourmash
----------------------------------------------------------------------------------------
*/

include { KRAKEN2              } from '../../modules/kraken2/main'
include { BRACKEN              } from '../../modules/bracken/main'
include { KRAKEN_BIOM          } from '../../modules/kraken_biom/main'
include { KRAKEN_TO_PHYLOSEQ   } from '../../modules/kraken_to_phyloseq/main'
include { TAX_REPORT_KRAKEN2   } from '../../modules/tax_report_kraken2/main'
include { SOURMASH             } from '../../modules/sourmash/main'
include { SOURMASH_TO_PHYLOSEQ } from '../../modules/sourmash_to_phyloseq/main'
include { TAX_REPORT_SOURMASH  } from '../../modules/tax_report_sourmash/main'

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
        ch_k2_taxonomy = KRAKEN2(
            reads.combine(kraken_db),
            params.kraken_db_used
        )

        // Generate taxonomy report
        TAX_REPORT_KRAKEN2(
            ch_k2_taxonomy.report
                .concat(reads_report)
                .collect()
        )

        // Bracken abundance estimation
        ch_taxonomy_estimation = BRACKEN(
            ch_k2_taxonomy.kraken.combine(kraken_db),
            params.kraken_db_used
        )

        // Generate BIOM file
        ch_kraken_biom = KRAKEN_BIOM(
            ch_taxonomy_estimation.collect(),
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
