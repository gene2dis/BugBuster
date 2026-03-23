/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE DATABASES SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Handles database downloads and formatting for all pipeline components
    Updated: 2026-01-29 - Single-pass decontamination optimization
----------------------------------------------------------------------------------------
*/

include { FORMAT_KRAKEN_DB        } from '../../modules/local/format_db/main'
include { FORMAT_BOWTIE_INDEX     } from '../../modules/local/format_db/main'
include { FORMAT_NT_BLAST_DB      } from '../../modules/local/format_db/main'
include { FORMAT_TAXDUMP_FILES    } from '../../modules/local/format_db/main'
include { DOWNLOAD_DEEPARG_DB     } from '../../modules/local/format_db/main'
include { FORMAT_CHECKM2_DB       } from '../../modules/local/format_db/main'
include { BUILD_PHIX_BOWTIE2_INDEX } from '../../modules/local/format_db/main'
include { DOWNLOAD_GTDBTK_DB      } from '../../modules/local/format_db/main'
include { SOURMASH_TAX_PREPARE    } from '../../modules/local/format_db/main'
include { RGI_LOAD                } from '../../modules/local/rgi_load/main'
include { RGI_LOAD_WILDCARD       } from '../../modules/local/rgi_load_wildcard/main'
include { BOWTIE2_BUILD_COMBINED  } from '../../modules/local/bowtie2_build_combined/main'

workflow PREPARE_DATABASES {
    
    main:
    // Initialize empty channels
    ch_kraken_db        = Channel.empty()
    ch_sourmash_db      = Channel.empty()
    ch_karga_db         = Channel.empty()
    ch_kargva_db        = Channel.empty()
    ch_deeparg_db       = Channel.empty()
    ch_blast_db         = Channel.empty()
    ch_taxdump          = Channel.empty()
    ch_gtdbtk_db        = Channel.empty()
    ch_checkm2_db       = Channel.empty()
    ch_rgi_card_db      = Channel.empty()

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
            // Custom database: expect list with [kmer_db, lineages_file]
            ch_sourmash_files = Channel.fromList(params.custom_sourmash_db)
                .map { filepath -> file(filepath, checkIfExists: true) }
                .collect()
            
            ch_sourmash_kmer = ch_sourmash_files.map { files -> files[0] }
            ch_sourmash_lineages = ch_sourmash_files.map { files -> files[1] }
        } else {
            // Reference database: download both k-mer and lineages files
            ch_sourmash_files = Channel.fromList(params.sourmash_ref_db[params.sourmash_db]["file"])
                .map { filepath -> file(filepath) }
                .collect()
            
            ch_sourmash_kmer = ch_sourmash_files.map { files -> files[0] }
            ch_sourmash_lineages = ch_sourmash_files.map { files -> files[1] }
        }
        
        // Prepare taxonomy database ONCE (not per-sample)
        ch_sourmash_tax_db = SOURMASH_TAX_PREPARE(ch_sourmash_lineages)
        
        // Combine k-mer DB and prepared taxonomy DB for downstream use
        ch_sourmash_db = ch_sourmash_kmer
            .combine(ch_sourmash_tax_db)
            .collect()
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
    // Combined decontamination index (PhiX + host) for QC
    //
    if ( params.quality_control ) {
        if ( params.custom_decontamination_index ) {
            // Use pre-built combined index
            ch_decontamination_index = Channel.fromPath(params.custom_decontamination_index, checkIfExists: true)
        } else {
            // Collect FASTA file paths into lists
            def phix_files = params.custom_phiX_fasta ? 
                [params.custom_phiX_fasta] : 
                params.bowtie_ref_genomes_for_build[params.phiX_index]["file"]
            
            def host_files = params.custom_host_fasta ? 
                [params.custom_host_fasta] : 
                params.bowtie_ref_host_index[params.host_db]["file"]
            
            // Combine lists and create single channel
            def all_fasta_files = phix_files + host_files
            
            // Build combined index from all FASTA files
            BOWTIE2_BUILD_COMBINED(
                Channel.fromList(all_fasta_files).map { filepath -> file(filepath) }.collect(),
                "contaminants"
            )
            
            // Extract only the index output (not versions)
            ch_decontamination_index = BOWTIE2_BUILD_COMBINED.out.index
        }
    } else {
        ch_decontamination_index = Channel.empty()
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

    //
    // RGI CARD database for AMR prediction
    //
    if ( params.rgi_prediction ) {
        if ( params.custom_rgi_card_db && params.custom_rgi_wildcard ) {
            // Use existing CARD database and add custom WildCARD
            ch_card_base = Channel.fromPath(params.custom_rgi_card_db, checkIfExists: true)
            ch_wildcard = Channel.fromPath(params.custom_rgi_wildcard, checkIfExists: true)
            ch_rgi_card_db = RGI_LOAD_WILDCARD(ch_card_base, ch_wildcard).card_db
        } else if ( params.custom_rgi_card_db ) {
            // Use existing pre-prepared CARD database (may or may not include WildCARD)
            ch_rgi_card_db = Channel.fromPath(params.custom_rgi_card_db, checkIfExists: true)
        } else {
            // Download and prepare CARD database (with optional WildCARD)
            ch_rgi_card_db = RGI_LOAD(
                params.rgi_card_version,
                params.rgi_include_wildcard
            ).card_db
        }
    }

    emit:
    kraken_db              = ch_kraken_db
    sourmash_db            = ch_sourmash_db
    karga_db               = ch_karga_db
    kargva_db              = ch_kargva_db
    decontamination_index  = ch_decontamination_index  // Combined phiX + host index
    deeparg_db             = ch_deeparg_db
    blast_db               = ch_blast_db
    taxdump                = ch_taxdump
    gtdbtk_db              = ch_gtdbtk_db
    checkm2_db             = ch_checkm2_db
    rgi_card_db            = ch_rgi_card_db
}
