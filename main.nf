#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gene2dis/BugBuster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/gene2dis/BugBuster
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

import groovy.transform.Field

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE LOGO AND INFO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

@Field
def logo = '''
\u001B[0m
     \u001B[31m╔███████████╗       \u001B[36m██████╗ ██╗   ██╗ ██████╗
   \u001B[31m╔██╝   \u001B[32m▄ ▄   \u001B[31m╚▀█╗\u001B[36m     ██╔══██╗██║   ██║██╔════╝
 \u001B[31m╔█▀▀╚▀█╗\u001B[32m▄▄█▄▄    \u001B[31m▀▀█╗\u001B[36m   ██████╔╝██║   ██║██║  ███╗
\u001B[31m██╝ \u001B[32m▄  ▄\u001B[31m█╗\u001B[33mo  o\u001B[32m▀▄  ▄ \u001B[31m╚██\u001B[36m  ██╔══██╗██║   ██║██║   ██║
\u001B[31m██   \u001B[32m▀▀█\u001B[31m╚▀█╗\u001B[33mo  \u001B[32m█▀▀   \u001B[31m██\u001B[36m  ██████╔╝╚██████╔╝╚██████╔╝
\u001B[31m██  \u001B[32m▄  █  \u001B[31m╚█▄╗ \u001B[32m█  ▄  \u001B[31m██\u001B[36m  ╚═════╝  ╚═════╝  ╚═════╝
\u001B[31m██  \u001B[32m▄▀▀█\u001B[33m o  \u001B[31m╚█▄\u001B[32m█▀▀▄\u001B[31m  ██\u001B[36m  ██████╗ ██╗   ██╗███████╗████████╗███████╗██████╗
\u001B[31m██     \u001B[32m█\u001B[33m   o  \u001B[31m╚█╗    ██\u001B[36m  ██╔══██╗██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗
\u001B[31m██╗  \u001B[32m▄▀▀▄\u001B[33m o  o\u001B[32m▄▀\u001B[31m▀█╗ ╔██\u001B[36m  ██████╔╝██║   ██║███████╗   ██║   █████╗  ██████╔╝
 \u001B[31m╚█▄▄    \u001B[32m▀▀█▀▀\u001B[0m   \u001B[31m╚█▄█╝\u001B[36m   ██╔══██╗██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗
   \u001B[31m╚█▄    \u001B[32m▀ ▀\u001B[0m    \u001B[31m▄█╝\u001B[36m     ██████╔╝╚██████╔╝███████║   ██║   ███████╗██║  ██║
     \u001B[31m╚███████████╝\u001B[36m       ╚═════╝  ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝
\u001B[0m
'''

def printVersion() {
    log.info ""
    log.info "  ${workflow.manifest.name} v${workflow.manifest.version}"
    log.info "  ${workflow.manifest.description}"
    log.info ""
}

def printHelp() {
    log.info logo
    log.info """
    \u001B[1;33mUsage:\u001B[0m

    The typical command for running the pipeline is as follows:

      nextflow run main.nf --input samplesheet.csv --output ./results -profile docker

    \u001B[1;33mMandatory arguments:\u001B[0m
      --input                       Path to CSV samplesheet with columns: sample,r1,r2,s
      --output                      Path to output directory

    \u001B[1;33mPipeline options:\u001B[0m
      --quality_control             Enable QC and host filtering (default: ${params.quality_control})
      --assembly_mode               Assembly mode: 'assembly', 'coassembly', 'none' (default: ${params.assembly_mode})
      --taxonomic_profiler          Profiler: 'kraken2', 'sourmash', 'none' (default: ${params.taxonomic_profiler})
      --include_binning             Enable binning and refinement (default: ${params.include_binning})
      --read_arg_prediction         Enable read-level ARG prediction (default: ${params.read_arg_prediction})
      --contig_tax_and_arg          Enable contig-level taxonomy and ARG (default: ${params.contig_tax_and_arg})
      --contig_level_metacerberus   Enable MetaCerberus annotation (default: ${params.contig_level_metacerberus})

    \u001B[1;33mResource options:\u001B[0m
      --max_cpus                    Maximum CPUs per process (default: ${params.max_cpus})
      --max_memory                  Maximum memory per process (default: ${params.max_memory})
      --max_time                    Maximum time per process (default: ${params.max_time})

    \u001B[1;33mProfile options:\u001B[0m
      -profile docker               Run with Docker containers
      -profile singularity          Run with Singularity containers
      -profile podman               Run with Podman containers
      -profile conda                Run with Conda environments
      -profile slurm_singularity    Run on SLURM with Singularity
      -profile test                 Run with minimal test dataset

    \u001B[1;33mOther options:\u001B[0m
      --help                        Show this help message
      --version                     Show pipeline version

    For more information, visit: ${workflow.manifest.homePage}
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS AND PRINT INFO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Show help message
if (params.help) {
    printHelp()
    exit 0
}

// Show version
if (params.containsKey('version') && params.version) {
    printVersion()
    exit 0
}

// Print logo
log.info logo
log.info ""
log.info "  ${workflow.manifest.name} v${workflow.manifest.version}"
log.info "  ================================================"
log.info ""

// Validate required parameters
if (!params.input) {
    log.error "ERROR: --input parameter is required"
    printHelp()
    exit 1
}

if (!params.output) {
    log.error "ERROR: --output parameter is required"
    printHelp()
    exit 1
}

// Print run configuration
log.info "  Run configuration:"
log.info "  -------------------"
log.info "  Input samplesheet    : ${params.input}"
log.info "  Output directory     : ${params.output}"
log.info "  Quality control      : ${params.quality_control}"
log.info "  Assembly mode        : ${params.assembly_mode}"
log.info "  Taxonomic profiler   : ${params.taxonomic_profiler}"
log.info "  Include binning      : ${params.include_binning}"
log.info "  Read ARG prediction  : ${params.read_arg_prediction}"
log.info "  Contig tax and ARG   : ${params.contig_tax_and_arg}"
log.info ""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Subworkflows
include { INPUT_CHECK        } from './subworkflows/local/input_check'
include { PREPARE_DATABASES  } from './subworkflows/local/prepare_databases'
include { QC                 } from './subworkflows/local/qc'
include { TAXONOMY           } from './subworkflows/local/taxonomy'
include { ASSEMBLY           } from './subworkflows/local/assembly'
include { BINNING            } from './subworkflows/local/binning'

// Modules for functionality not covered by subworkflows

	// FUNCTIONAL ANNOTATION
include { METACERBERUS_CONTIGS } from './modules/metacerberus/main'

	// TAXONOMIC PREDICTION IN CONTIGS
include { NT_BLASTN        } from './modules/nt_blastn/main'
include { BLOBTOOLS        } from './modules/blobtools/main'
include { SAMTOOLS_INDEX   } from './modules/samtools_index/main'
include { BLOBPLOT         } from './modules/blobplot/main'

	// ORF PREDICTION IN CONTIGS AND BINS
include { PRODIGAL_BINS    } from './modules/prodigal/main'
include { PRODIGAL_CONTIGS } from './modules/prodigal/main'

	// ARG PREDICTION IN READS
include { KARGVA           } from './modules/kargva/main'
include { KARGA            } from './modules/karga/main'
include { ARGS_OAP         } from './modules/args_oap/main'
include { ARG_NORM_REPORT  } from './modules/arg_norm_report/main'

	// ARG PREDICTION IN CONTIGS AND BINS
include { DEEPARG_BINS             } from './modules/deeparg/main'
include { DEEPARG_CONTIGS          } from './modules/deeparg/main'
include { ARG_CONTIG_LEVEL_REPORT  } from './modules/arg_contig_level_report/main'
include { ARG_FASTA_FORMATTER      } from './modules/arg_fasta_formatter/main'
include { CLUSTERING               } from './modules/clustering/main'
include { ARG_BLOBPLOT             } from './modules/arg_blobplot/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    //
    // SUBWORKFLOW: Validate and parse input samplesheet
    //
    INPUT_CHECK(file(params.input, checkIfExists: true))
    ch_reads = INPUT_CHECK.out.reads

    //
    // SUBWORKFLOW: Prepare all databases
    //
    PREPARE_DATABASES()
    
    //
    // SUBWORKFLOW: Quality control and host decontamination
    //
    QC(
        ch_reads,
        PREPARE_DATABASES.out.phix_index,
        PREPARE_DATABASES.out.host_index
    )
    
    ch_clean_reads           = QC.out.reads
    ch_clean_reads_coassembly = QC.out.reads_coassembly
    ch_reads_report          = QC.out.report

    //
    // SUBWORKFLOW: Taxonomic profiling
    //
    if ( params.taxonomic_profiler != "none" ) {
        TAXONOMY(
            ch_clean_reads,
            ch_reads_report,
            PREPARE_DATABASES.out.kraken_db,
            PREPARE_DATABASES.out.sourmash_db
        )
    }

    //
    // ARG PREDICTION IN READS
    //
    if ( params.read_arg_prediction ) {
        ch_args_oap = ARGS_OAP(ch_clean_reads)
        ch_argv_prediction = KARGVA(ch_clean_reads.combine(PREPARE_DATABASES.out.kargva_db))
        ch_arg_prediction = KARGA(ch_argv_prediction.kargva_reads.combine(PREPARE_DATABASES.out.karga_db))
        ARG_NORM_REPORT(
            ch_arg_prediction
                .concat(ch_argv_prediction.kargva_reports)
                .concat(ch_args_oap)
                .collect()
        )
    }

    //
    // SUBWORKFLOW: Assembly
    //
    ch_contigs_meta = Channel.empty()
    ch_bam_meta     = Channel.empty()
    ch_refined_bins = Channel.empty()
    
    if ( params.assembly_mode != "none" ) {
        ASSEMBLY(
            ch_clean_reads,
            ch_clean_reads_coassembly
        )
        
        ch_contigs_meta = ASSEMBLY.out.contigs_meta
        ch_bam_meta     = ASSEMBLY.out.bam_meta

        //
        // SUBWORKFLOW: Binning
        //
        if ( params.include_binning ) {
            BINNING(
                ch_clean_reads,
                ASSEMBLY.out.contigs,
                ASSEMBLY.out.bam,
                PREPARE_DATABASES.out.checkm2_db,
                PREPARE_DATABASES.out.gtdbtk_db
            )
            ch_refined_bins = BINNING.out.refined_bins
        }

        //
        // MetaCerberus annotation (per-sample assembly only)
        //
        if ( params.assembly_mode == "assembly" && params.contig_level_metacerberus ) {
            METACERBERUS_CONTIGS(ch_contigs_meta)
        }
    }

    //
    // CONTIG-LEVEL TAXONOMY AND ARG PREDICTION
    //
    if ( params.contig_tax_and_arg && params.assembly_mode != "none" ) {
        ch_nt_blastn = NT_BLASTN(ch_contigs_meta.combine(PREPARE_DATABASES.out.blast_db))
        ch_index_bam = SAMTOOLS_INDEX(ch_bam_meta)
        ch_blob_table = BLOBTOOLS(
            ch_nt_blastn
                .join(ch_index_bam)
                .combine(PREPARE_DATABASES.out.taxdump.collect())
        )
        BLOBPLOT(ch_blob_table.only_blob.collect())

        ch_contig_proteins = PRODIGAL_CONTIGS(ch_contigs_meta)
        ch_contig_args = DEEPARG_CONTIGS(ch_contig_proteins.combine(PREPARE_DATABASES.out.deeparg_db))
        ch_arg_contig_data = ARG_CONTIG_LEVEL_REPORT(
            ch_contig_args.only_deeparg
                .concat(ch_blob_table.only_blob)
                .collect()
        )
        ARG_BLOBPLOT(ch_arg_contig_data)
    }

    //
    // ARG PREDICTION IN BINS AND CLUSTERING
    //
    if ( params.arg_bin_clustering && params.include_binning ) {
        ch_raw_orfs = PRODIGAL_BINS(ch_refined_bins)
        ch_deeparg = DEEPARG_BINS(ch_raw_orfs.combine(PREPARE_DATABASES.out.deeparg_db))
        ch_arg_fasta = ARG_FASTA_FORMATTER(ch_raw_orfs.join(ch_deeparg))
        ch_clusters = CLUSTERING(ch_arg_fasta.collect())
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        Exit status  : ${workflow.exitStatus}
        Work dir     : ${workflow.workDir}
        Output dir   : ${params.output}
        """
        .stripIndent()

    log.info msg

    if (workflow.success) {
        log.info "\u001B[32m========================================\u001B[0m"
        log.info "\u001B[32m  Pipeline completed successfully!\u001B[0m"
        log.info "\u001B[32m========================================\u001B[0m"
    } else {
        log.error "\u001B[31m========================================\u001B[0m"
        log.error "\u001B[31m  Pipeline completed with errors\u001B[0m"
        log.error "\u001B[31m========================================\u001B[0m"
    }
}

workflow.onError {
    log.error "Pipeline failed. Check error message above or in ${params.output}/pipeline_info/"
}
