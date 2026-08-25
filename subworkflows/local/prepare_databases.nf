/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE DATABASES SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Handles database downloads and formatting for all pipeline components
    Updated: 2026-01-29 - Single-pass decontamination optimization
----------------------------------------------------------------------------------------
*/

include { FORMAT_KRAKEN_DB        } from '../../modules/local/format_db/main'
include { FORMAT_NT_BLAST_DB      } from '../../modules/local/format_db/main'
include { FORMAT_TAXDUMP_FILES    } from '../../modules/local/format_db/main'
include { DOWNLOAD_DEEPARG_DB     } from '../../modules/local/format_db/main'
include { FORMAT_CHECKM2_DB       } from '../../modules/local/format_db/main'
include { DOWNLOAD_GTDBTK_DB      } from '../../modules/local/format_db/main'
include { SOURMASH_TAX_PREPARE    } from '../../modules/local/format_db/main'
include { RGI_LOAD                } from '../../modules/local/rgi_load/main'
include { RGI_LOAD_WILDCARD       } from '../../modules/local/rgi_load_wildcard/main'
include { BOWTIE2_BUILD_COMBINED  } from '../../modules/local/bowtie2_build_combined/main'

workflow PREPARE_DATABASES {
    
    main:
    // Initialize empty channels
    ch_kraken_db        = channel.empty()
    ch_sourmash_db      = channel.empty()
    ch_karga_db         = channel.empty()
    ch_kargva_db        = channel.empty()
    ch_deeparg_db       = channel.empty()
    ch_blast_db         = channel.empty()
    ch_taxdump          = channel.empty()
    ch_gtdbtk_db        = channel.empty()
    ch_checkm2_db       = channel.empty()
    ch_rgi_card_db      = channel.empty()
    ch_versions         = channel.empty()

    //
    // Kraken2 database
    //
    if ( params.taxonomic_profiler == "kraken2" ) {
        if ( params.custom_kraken_db ) {
            ch_kraken_db = channel.fromPath(params.custom_kraken_db, checkIfExists: true)
        } else {
            ch_kraken_ref = channel.fromList(params.kraken_ref_db[params.kraken2_db]["url"])
            ch_kraken_db = FORMAT_KRAKEN_DB(ch_kraken_ref)
        }
    }

    //
    // Sourmash database
    //
    if ( params.taxonomic_profiler == "sourmash" ) {
        if ( params.custom_sourmash_db ) {
            // Custom database: expect list with [kmer_db, lineages_file]
            ch_sourmash_files = channel.fromList(params.custom_sourmash_db)
                .map { filepath -> file(filepath, checkIfExists: true) }
                .collect()
            
            ch_sourmash_kmer = ch_sourmash_files.map { files -> files[0] }
            ch_sourmash_lineages = ch_sourmash_files.map { files -> files[1] }
        } else {
            // Reference database: download both k-mer and lineages files
            ch_sourmash_files = channel.fromList(params.sourmash_ref_db[params.sourmash_db]["url"])
                .map { filepath -> file(filepath) }
                .collect()
            
            ch_sourmash_kmer = ch_sourmash_files.map { files -> files[0] }
            ch_sourmash_lineages = ch_sourmash_files.map { files -> files[1] }
        }
        
        // Prepare taxonomy database ONCE (not per-sample)
        SOURMASH_TAX_PREPARE(ch_sourmash_lineages)
        ch_sourmash_tax_db = SOURMASH_TAX_PREPARE.out.tax_db
        ch_versions = ch_versions.mix(SOURMASH_TAX_PREPARE.out.versions)
        
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
            ch_karga_db = channel.of(file(params.custom_karga_db, checkIfExists: true))
        } else {
            ch_karga_db = channel.fromList(params.karga_ref_db[params.karga_db]["url"])
                .map { filepath -> file(filepath) }
        }

        if ( params.custom_kargva_db ) {
            ch_kargva_db = channel.of(file(params.custom_kargva_db, checkIfExists: true))
        } else {
            ch_kargva_db = channel.fromList(params.kargva_ref_db[params.kargva_db]["url"])
                .map { filepath -> file(filepath) }
        }
    }

    //
    // Combined decontamination index (PhiX + host) for QC
    //
    if ( params.quality_control ) {
        if ( params.custom_decontamination_index ) {
            // Use pre-built combined index
            ch_decontamination_index = channel.fromPath(params.custom_decontamination_index, checkIfExists: true)
        } else {
            // Collect FASTA file paths into lists
            def phix_files = params.custom_phiX_fasta ? 
                [params.custom_phiX_fasta] : 
                params.bowtie_ref_genomes_for_build[params.phiX_index]["url"]
            
            def host_files = params.custom_host_fasta ? 
                [params.custom_host_fasta] : 
                params.bowtie_ref_host_index[params.host_db]["url"]
            
            // Combine lists and create single channel
            def all_fasta_files = phix_files + host_files

            // These entries must be genome FASTA (optionally gzipped) — a prebuilt
            // bowtie2 index cannot be concatenated and rebuilt
            def non_fasta = all_fasta_files.findAll { f ->
                f.toString() ==~ /.*\.(zip|bt2l?|tar|tar\.gz|tgz)$/
            }
            if ( non_fasta ) {
                error "Decontamination references must be genome FASTA files, but got: ${non_fasta.join(', ')}. " +
                      "To use a pre-built Bowtie2 index, pass it via --custom_decontamination_index instead."
            }

            // Build combined index from all FASTA files
            BOWTIE2_BUILD_COMBINED(
                channel.fromList(all_fasta_files).map { filepath -> file(filepath) }.collect(),
                "contaminants"
            )
            
            // Extract only the index output (not versions)
            ch_decontamination_index = BOWTIE2_BUILD_COMBINED.out.index
            ch_versions = ch_versions.mix(BOWTIE2_BUILD_COMBINED.out.versions)
        }
    } else {
        ch_decontamination_index = channel.empty()
    }

    //
    // DeepARG for contig-level analysis and/or bin-level ARG clustering
    //
    if ( params.contig_tax_and_arg || params.arg_bin_clustering ) {
        if ( params.custom_deeparg_db ) {
            ch_deeparg_db = channel.fromPath(params.custom_deeparg_db, checkIfExists: true)
        } else {
            DOWNLOAD_DEEPARG_DB()
            ch_deeparg_db = DOWNLOAD_DEEPARG_DB.out.deeparg_db
            ch_versions = ch_versions.mix(DOWNLOAD_DEEPARG_DB.out.versions)
        }
    }

    //
    // BLAST and taxdump for contig-level taxonomy
    //
    if ( params.contig_tax_and_arg ) {
        if ( params.custom_blast_db ) {
            ch_blast_db = channel.fromPath(params.custom_blast_db, checkIfExists: true)
        } else {
            ch_blast_ref = channel.fromList(params.blast_ref_db[params.blast_db]["url"])
            ch_blast_db = FORMAT_NT_BLAST_DB(ch_blast_ref)
        }

        if ( params.custom_taxdump_files ) {
            ch_taxdump = channel.fromPath(params.custom_taxdump_files, checkIfExists: true)
        } else {
            ch_taxdump_ref = channel.fromList(params.taxonomy_files[params.taxdump_files]["url"])
            ch_taxdump = FORMAT_TAXDUMP_FILES(ch_taxdump_ref)
        }
    }

    //
    // GTDB-TK and CheckM2 databases for binning
    //
    if ( params.include_binning ) {
        if ( params.custom_gtdbtk_db ) {
            ch_gtdbtk_db = channel.fromPath(params.custom_gtdbtk_db, checkIfExists: true)
        } else {
            ch_gtdbtk_ref = channel.fromList(params.gtdbtk_ref_db[params.gtdbtk_db]["url"])
            ch_gtdbtk_db = DOWNLOAD_GTDBTK_DB(ch_gtdbtk_ref)
        }

        if ( params.custom_checkm2_db ) {
            ch_checkm2_db = channel.fromPath(params.custom_checkm2_db, checkIfExists: true)
        } else {
            ch_checkm2_ref = channel.fromList(params.checkm2_ref_db[params.checkm2_db]["url"])
            ch_checkm2_db = FORMAT_CHECKM2_DB(ch_checkm2_ref)
        }
    }

    //
    // RGI CARD database for AMR prediction
    //
    if ( params.rgi_prediction ) {
        if ( params.custom_rgi_card_db && params.custom_rgi_wildcard ) {
            // Use existing CARD database and add custom WildCARD
            ch_card_base = channel.fromPath(params.custom_rgi_card_db, checkIfExists: true)
            ch_wildcard = channel.fromPath(params.custom_rgi_wildcard, checkIfExists: true)
            ch_rgi_card_db = RGI_LOAD_WILDCARD(ch_card_base, ch_wildcard).card_db
            ch_versions = ch_versions.mix(RGI_LOAD_WILDCARD.out.versions)
        } else if ( params.custom_rgi_card_db ) {
            // Use existing pre-prepared CARD database (may or may not include WildCARD)
            ch_rgi_card_db = channel.fromPath(params.custom_rgi_card_db, checkIfExists: true)
        } else {
            // Download and prepare CARD database (with optional WildCARD)
            ch_rgi_card_db = RGI_LOAD(
                params.rgi_card_version,
                params.rgi_include_wildcard
            ).card_db
            ch_versions = ch_versions.mix(RGI_LOAD.out.versions)
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
    versions               = ch_versions
}
