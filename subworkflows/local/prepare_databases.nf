/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE DATABASES SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Handles database downloads and formatting for all pipeline components
----------------------------------------------------------------------------------------
*/

include { FORMAT_KRAKEN_DB        } from '../../modules/format_db/main'
include { FORMAT_BOWTIE_INDEX     } from '../../modules/format_db/main'
include { FORMAT_NT_BLAST_DB      } from '../../modules/format_db/main'
include { FORMAT_TAXDUMP_FILES    } from '../../modules/format_db/main'
include { DOWNLOAD_DEEPARG_DB     } from '../../modules/format_db/main'
include { FORMAT_CHECKM2_DB       } from '../../modules/format_db/main'
include { BUILD_PHIX_BOWTIE2_INDEX } from '../../modules/format_db/main'
include { DOWNLOAD_GTDBTK_DB      } from '../../modules/format_db/main'

workflow PREPARE_DATABASES {
    
    main:
    // Initialize empty channels
    ch_kraken_db        = Channel.empty()
    ch_sourmash_db      = Channel.empty()
    ch_karga_db         = Channel.empty()
    ch_kargva_db        = Channel.empty()
    ch_phix_index       = Channel.empty()
    ch_host_index       = Channel.empty()
    ch_deeparg_db       = Channel.empty()
    ch_blast_db         = Channel.empty()
    ch_taxdump          = Channel.empty()
    ch_gtdbtk_db        = Channel.empty()
    ch_checkm2_db       = Channel.empty()

    //
    // Kraken2 database
    //
    if ( params.taxonomic_profiler == "kraken2" ) {
        if ( params.custom_kraken_db ) {
            ch_kraken_db = Channel.fromPath(params.custom_kraken_db, checkIfExists: true)
        } else {
            ch_kraken_ref = Channel.fromList(params.kraken_ref_db[params.kraken2_db]["file"])
            ch_kraken_db = FORMAT_KRAKEN_DB(ch_kraken_ref)
        }
    }

    //
    // Sourmash database
    //
    if ( params.taxonomic_profiler == "sourmash" ) {
        if ( params.custom_sourmash_db ) {
            ch_sourmash_db = Channel.fromList(params.custom_sourmash_db)
                .map { filepath -> file(filepath, checkIfExists: true) }
                .collect()
        } else {
            ch_sourmash_db = Channel.fromList(params.sourmash_ref_db[params.sourmash_db]["file"])
                .map { filepath -> file(filepath) }
        }
    }

    //
    // KARGA/KARGVA databases for read-level ARG prediction
    //
    if ( params.read_arg_prediction ) {
        if ( params.custom_karga_db ) {
            ch_karga_db = Channel.of(file(params.custom_karga_db, checkIfExists: true))
        } else {
            ch_karga_db = Channel.fromList(params.karga_ref_db[params.karga_db]["file"])
                .map { filepath -> file(filepath) }
        }

        if ( params.custom_kargva_db ) {
            ch_kargva_db = Channel.of(file(params.custom_kargva_db, checkIfExists: true))
        } else {
            ch_kargva_db = Channel.fromList(params.kargva_ref_db[params.kargva_db]["file"])
                .map { filepath -> file(filepath) }
        }
    }

    //
    // PhiX and host Bowtie2 indexes for QC
    //
    if ( params.quality_control ) {
        if ( params.custom_phiX_index ) {
            ch_phix_index = Channel.fromPath(params.custom_phiX_index, checkIfExists: true)
        } else {
            ch_phix_ref = Channel.fromList(params.bowtie_ref_genomes_for_build[params.phiX_index]["file"])
                .map { filepath -> file(filepath) }
            ch_phix_index = BUILD_PHIX_BOWTIE2_INDEX(ch_phix_ref)
        }

        if ( params.custom_bowtie_host_index ) {
            ch_host_index = Channel.fromPath(params.custom_bowtie_host_index, checkIfExists: true)
        } else {
            ch_host_ref = Channel.fromList(params.bowtie_ref_host_index[params.host_db]["file"])
            ch_host_index = FORMAT_BOWTIE_INDEX(ch_host_ref)
        }
    }

    //
    // DeepARG, BLAST, and taxdump for contig-level analysis
    //
    if ( params.contig_tax_and_arg ) {
        if ( params.custom_deeparg_db ) {
            ch_deeparg_db = Channel.fromPath(params.custom_deeparg_db, checkIfExists: true)
        } else {
            ch_deeparg_db = DOWNLOAD_DEEPARG_DB()
        }

        if ( params.custom_blast_db ) {
            ch_blast_db = Channel.fromPath(params.custom_blast_db, checkIfExists: true)
        } else {
            ch_blast_ref = Channel.fromList(params.blast_ref_db[params.blast_db]["url"])
            ch_blast_db = FORMAT_NT_BLAST_DB(ch_blast_ref)
        }

        if ( params.custom_taxdump_files ) {
            ch_taxdump = Channel.fromPath(params.custom_taxdump_files, checkIfExists: true)
        } else {
            ch_taxdump_ref = Channel.fromList(params.taxonomy_files[params.taxdump_files]["url"])
            ch_taxdump = FORMAT_TAXDUMP_FILES(ch_taxdump_ref)
        }
    }

    //
    // GTDB-TK and CheckM2 databases for binning
    //
    if ( params.include_binning ) {
        if ( params.custom_gtdbtk_db ) {
            ch_gtdbtk_db = Channel.fromPath(params.custom_gtdbtk_db, checkIfExists: true)
        } else {
            ch_gtdbtk_ref = Channel.fromList(params.gtdbtk_ref_db[params.gtdbtk_db]["url"])
            ch_gtdbtk_db = DOWNLOAD_GTDBTK_DB(ch_gtdbtk_ref)
        }

        if ( params.custom_checkm2_db ) {
            ch_checkm2_db = Channel.fromPath(params.custom_checkm2_db, checkIfExists: true)
        } else {
            ch_checkm2_ref = Channel.fromList(params.checkm2_ref_db[params.checkm2_db]["url"])
            ch_checkm2_db = FORMAT_CHECKM2_DB(ch_checkm2_ref)
        }
    }

    emit:
    kraken_db   = ch_kraken_db
    sourmash_db = ch_sourmash_db
    karga_db    = ch_karga_db
    kargva_db   = ch_kargva_db
    phix_index  = ch_phix_index
    host_index  = ch_host_index
    deeparg_db  = ch_deeparg_db
    blast_db    = ch_blast_db
    taxdump     = ch_taxdump
    gtdbtk_db   = ch_gtdbtk_db
    checkm2_db  = ch_checkm2_db
}
