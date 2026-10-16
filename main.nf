#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gene2dis/BugBuster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/gene2dis/BugBuster
----------------------------------------------------------------------------------------
*/

include { validateParameters } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE LOGO AND INFO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def pipelineLogo() {
    return '''
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
}

def printVersion() {
    log.info ""
    log.info "  ${workflow.manifest.name} v${workflow.manifest.version}"
    log.info "  ${workflow.manifest.description}"
    log.info ""
}

def printHelp() {
    log.info pipelineLogo()
    log.info """
    \u001B[1;33mUsage:\u001B[0m

    The typical command for running the pipeline is as follows:

      nextflow run main.nf --input samplesheet.csv --output ./results -profile docker

    \u001B[1;33mMandatory arguments:\u001B[0m
      --input                       Path to CSV samplesheet with columns: sample,r1,r2,s
      --output                      Path to output directory

    \u001B[1;33mPipeline options:\u001B[0m
      --quality_control             Enable QC and host filtering (default: ${params.quality_control})
      --read_aligner                Aligner for host removal, reads->contigs and bin depth: 'bwamem3' (needs AVX2), 'bowtie2' (default: ${params.read_aligner})
      --assembly_mode               Assembly mode: 'assembly', 'coassembly', 'none' (default: ${params.assembly_mode})
      --taxonomic_profiler          Profiler: 'kraken2', 'sourmash', 'sylph', 'none' (default: ${params.taxonomic_profiler})
      --include_binning             Enable binning and refinement (default: ${params.include_binning})
      --binners                     Comma-separated binners: comebin, semibin, metabat2 (default: ${params.binners})
      --read_arg_prediction         Enable read-level ARG prediction (default: ${params.read_arg_prediction})
      --rgi_prediction              Enable RGI AMR prediction with pathogen-of-origin (default: ${params.rgi_prediction})
      --kmarvel_prediction          Enable K-MARVEL read-level ARG assembly against CARD (default: ${params.kmarvel_prediction})
      --contig_tax_and_arg          Enable contig-level taxonomy and ARG (default: ${params.contig_tax_and_arg})
      --contig_level_functional     Enable contig-level functional annotation (eggNOG-mapper v3 + run_dbcan) (default: ${params.contig_level_functional})
      --microbecensus               Estimate average genome size for CPGE normalization (contig/read functional branches; default: ${params.microbecensus})
      --functional_cazy             Run run_dbcan CAZy annotation on predicted proteins (functional branch; default: ${params.functional_cazy})
      --dbcan_consensus             dbCAN calls feeding the summary tables: recommended | any (default: ${params.dbcan_consensus})
      --mag_level_functional        Bakta annotation of refined bins; needs --include_binning and >= 2 --binners (default: ${params.mag_level_functional})
      --read_level_functional       Read-level functional profiling backend: 'woltka', 'superfocus', 'humann', 'none' (default: ${params.read_level_functional})
      --woltka_uniq                 Woltka: leave multi-hit reads unassigned instead of dividing them 1/k (default: ${params.woltka_uniq})
      --superfocus_aligner          SUPER-FOCUS search backend: 'diamond', 'mmseqs2' (default: ${params.superfocus_aligner})

    \u001B[1;33mResults summary (00_summary/):\u001B[0m
      --results_summary             Tables, Excel workbook, figures and HTML report of QC, taxonomy and K-MARVEL results (default: ${params.results_summary})
      --sample_metadata             TSV/CSV with a 'sample' column and any metadata columns (default: none)
      --summary_group_by            Metadata column that groups samples, e.g. subject (default: none)
      --summary_order_by            Metadata column that orders samples within a group, e.g. visit (default: none)
      --summary_batch_label         Label in the summary's 'batch' column (default: the output directory's name)
      --summary_min_clean_pairs     Flag samples below this many clean read pairs as low depth (default: ${params.summary_min_clean_pairs})
      --card_prevalence             WildCARD card_prevalence.txt, adds candidate ARG hosts (default: none)

    \u001B[1;33mResource options:\u001B[0m
      --max_cpus                    Maximum CPUs per process (default: ${params.max_cpus})
      --max_memory                  Maximum memory per process (default: ${params.max_memory})
      --max_time                    Maximum time per process (default: ${params.max_time})

    \u001B[1;33mProfile options:\u001B[0m
      -profile docker               Run with Docker containers
      -profile singularity          Run with Singularity containers
      -profile podman               Run with Podman containers
      -profile apptainer            Run with Apptainer containers
      -profile conda                Run with Conda environments (containers recommended)
      -profile slurm                Run on a SLURM HPC cluster (combine: slurm,singularity)
      -profile aws                  Run on AWS Batch (combine: aws,docker)
      -profile gcp                  Run on Google Cloud Batch (combine: gcp,docker)
      -profile azure                Run on Azure Batch (combine: azure,docker)
      -profile low_disk             Delete the work dir after a successful run (completed runs not resumable)
      -profile test                 Run with minimal test dataset

    \u001B[1;33mOther options:\u001B[0m
      --help                        Show this help message
      --version                     Show pipeline version

    For more information, visit: ${workflow.manifest.homePage}
    """
}

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
include { FUNCTIONAL_ANNOTATION } from './subworkflows/local/functional_annotation'
include { READ_FUNCTIONAL    } from './subworkflows/local/read_functional'

// Modules for functionality not covered by subworkflows

	// FUNCTIONAL ANNOTATION
include { MICROBECENSUS } from './modules/local/microbecensus/main'

	// TAXONOMIC PREDICTION IN CONTIGS
include { NT_BLASTN        } from './modules/local/nt_blastn/main'
include { BLOBTOOLS        } from './modules/local/blobtools/main'
include { SAMTOOLS_INDEX as NFCORE_SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { BLOBPLOT         } from './modules/local/blobplot/main'

	// ORF PREDICTION IN CONTIGS AND BINS
	// PYRODIGAL is the shared contig gene-calling step (design doc Q9): one
	// pass feeds both DeepARG (contig_tax_and_arg) and functional annotation
	// (contig_level_functional)
include { PRODIGAL_BINS    } from './modules/local/prodigal/main'
include { PYRODIGAL        } from './modules/nf-core/pyrodigal/main'

	// ARG PREDICTION IN READS
include { KARGVA           } from './modules/local/kargva/main'
include { KARGA            } from './modules/local/karga/main'
include { ARGS_OAP         } from './modules/local/args_oap/main'
include { ARG_NORM_REPORT  } from './modules/local/arg_norm_report/main'

	// RGI AMR PREDICTION
include { RGI_BWT          } from './modules/local/rgi_bwt/main'
include { RGI_KMER         } from './modules/local/rgi_kmer/main'
include { RGI_REPORT       } from './modules/local/rgi_report/main'

	// K-MARVEL AMR PREDICTION
include { KMARVEL          } from './modules/local/kmarvel/main'
include { KMARVEL_REPORT   } from './modules/local/kmarvel_report/main'

	// RUN-LEVEL SUMMARY (00_summary/)
include { RESULTS_SUMMARY  } from './modules/local/results_summary/main'

	// ARG PREDICTION IN CONTIGS AND BINS
include { DEEPARG_BINS             } from './modules/local/deeparg/main'
include { DEEPARG_CONTIGS          } from './modules/local/deeparg/main'
include { ARG_CONTIG_LEVEL_REPORT  } from './modules/local/arg_contig_level_report/main'
include { ARG_FASTA_FORMATTER      } from './modules/local/arg_fasta_formatter/main'
include { CLUSTERING               } from './modules/local/clustering/main'
include { flagOn; sylphDbLabel     } from './subworkflows/local/utils_params'
include { ARG_BLOBPLOT             } from './modules/local/arg_blobplot/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    //
    // Help / version, parameter validation, and startup banner. Strict syntax
    // does not allow top-level statements, so this runs first in the workflow.
    //

    // Show help message
    if (flagOn(params.help)) {
        printHelp()
        System.exit(0)
    }

    // Show version
    if (params.containsKey('version') && params.version) {
        printVersion()
        System.exit(0)
    }

    // Validate all parameters against nextflow_schema.json (nf-schema plugin):
    // required params (--input/--output), types, enums (assembly_mode,
    // taxonomic_profiler, ...). Cross-parameter rules the schema cannot express
    // are checked by hand below.
    validateParameters()

    // Print logo
    log.info pipelineLogo()
    log.info ""
    log.info "  ${workflow.manifest.name} v${workflow.manifest.version}"
    log.info "  ================================================"
    log.info ""

    // Cloud profiles need an object-storage work directory; there is no sane
    // default, so fail fast instead of falling back to a local ./work that the
    // cloud executor cannot use (audit #24)
    def cloud_workdir_requirements = [
        aws  : ['aws_workdir', 's3://my-bucket/work', 's3://'],
        gcp  : ['gcp_workdir', 'gs://my-bucket/work', 'gs://'],
        azure: ['azure_workdir', 'az://my-container/work', 'az://'],
    ]
    def active_profiles = workflow.profile.tokenize(',')
    cloud_workdir_requirements.each { profile_name, req ->
        def (param_name, example, scheme) = req
        def workdir_ok = params[param_name] || workflow.workDir.toString().startsWith(scheme)
        if (active_profiles.contains(profile_name) && !workdir_ok) {
            error("-profile ${profile_name} requires an object-storage work directory: pass --${param_name} ${example} (or -work-dir ${example})")
        }
    }

    // Parse and validate binners parameter
    def binners_list = params.binners instanceof List ? params.binners : params.binners.toString().tokenize(',').collect { b -> b.trim().toLowerCase() }
    def valid_binners = ['comebin', 'semibin', 'metabat2']
    def invalid_binners = binners_list.findAll { b -> !(b in valid_binners) }
    if (binners_list.isEmpty()) {
        error("--binners must specify at least one binner. Valid options: ${valid_binners.join(', ')}")
    }
    if (invalid_binners) {
        error("Invalid binner(s): ${invalid_binners.join(', ')}. Valid options: ${valid_binners.join(', ')}")
    }

    // Reject contradictory feature combinations instead of silently skipping stages
    if (flagOn(params.include_binning) && params.assembly_mode == 'none') {
        error("--include_binning requires an assembly (--assembly_mode assembly or coassembly), but --assembly_mode is 'none'")
    }
    if (flagOn(params.contig_tax_and_arg) && params.assembly_mode == 'none') {
        error("--contig_tax_and_arg requires an assembly (--assembly_mode assembly or coassembly), but --assembly_mode is 'none'")
    }
    if (flagOn(params.contig_level_functional) && params.assembly_mode == 'none') {
        error("--contig_level_functional requires an assembly (--assembly_mode assembly or coassembly), but --assembly_mode is 'none'")
    }
    if (flagOn(params.arg_bin_clustering) && !flagOn(params.include_binning)) {
        error("--arg_bin_clustering requires --include_binning (it runs on the refined bins)")
    }
    if (flagOn(params.mag_level_functional) && !flagOn(params.include_binning)) {
        error("--mag_level_functional requires --include_binning (Bakta annotates the refined bins)")
    }
    if (flagOn(params.mag_level_functional) && binners_list.size() < 2) {
        error("--mag_level_functional requires at least two --binners: MetaWRAP refinement and its completeness/contamination quality filter only run with >= 2 binners, and Bakta must only annotate quality-filtered bins. Got: ${binners_list.join(', ')}")
    }
    // Read-level functional backend (design doc Section 7: unknown values
    // are rejected at launch). The schema enum already enforces this; kept
    // here so the accepted set lives next to the other launch validations
    def valid_read_functional = ['woltka', 'superfocus', 'humann', 'none']
    if (!(params.read_level_functional in valid_read_functional)) {
        error("Invalid --read_level_functional '${params.read_level_functional}'. Valid options: ${valid_read_functional.join(', ')}")
    }
    // The summary's grouping/ordering columns live in the metadata table
    if ((params.summary_group_by || params.summary_order_by) && !params.sample_metadata) {
        error("--summary_group_by and --summary_order_by name columns of --sample_metadata: pass a metadata table (TSV/CSV with a 'sample' column) as well")
    }
    // A sylph database and its sylph-tax metadata must come from the same
    // release: a custom database with the registry's metadata (or the reverse)
    // would get wrong or default lineages from sylph-tax without failing
    if (params.taxonomic_profiler == 'sylph' && (params.custom_sylph_db as boolean) != (params.custom_sylph_taxonomy as boolean)) {
        error("--custom_sylph_db and --custom_sylph_taxonomy must be given together: a sylph database and the sylph-tax metadata of the same release (e.g. gtdb-r232-c200-dbv2.syl2db with gtdb_r232_metadata.tsv.gz). Pass both, or neither to download the --sylph_db entry")
    }
    def valid_superfocus_aligners = ['diamond', 'mmseqs2']
    if (!(params.superfocus_aligner in valid_superfocus_aligners)) {
        error("Invalid --superfocus_aligner '${params.superfocus_aligner}'. Valid options: ${valid_superfocus_aligners.join(', ')}")
    }

    // low_disk deletes the work dir when the run succeeds, so a completed run
    // cannot be resumed or extended with -resume (an interrupted one can) — a bad
    // pairing with the long functional annotation runs and their large
    // databases, where -resume matters most (design doc Q10: warn, stay
    // results-neutral)
    if (workflow.profile.tokenize(',').contains('low_disk') && (flagOn(params.contig_level_functional) || flagOn(params.mag_level_functional) || params.read_level_functional != 'none')) {
        log.warn "functional annotation under -profile low_disk: the work dir is deleted when the run succeeds, so a completed run cannot be resumed (re-running with changed options repeats the long eggNOG-mapper / read-level steps), and the functional databases (eggNOG 7 ~44 GB, dbCAN ~7.4 GB, Bakta full ~31.9 GB download / light ~1.3 GB, WoLr2 ~94 GB, SUPER-FOCUS DB_90 ~0.7-0.9 GB download, HUMAnN ~71 GB with the full ChocoPhlAn / ~33 GB EC-filtered) are stored at --databases_dir regardless of this profile"
    }

    // bwa-mem3 (the default aligner) builds the human decontamination index
    // with ~57-75 GB RAM, bowtie2 with ~6 GB. Under a lower --max_memory the
    // build stops with exit 3: say so before the downloads start
    if (params.read_aligner == 'bwamem3' && flagOn(params.quality_control) && !params.custom_decontamination_index && !params.custom_host_fasta && (params.max_memory as MemoryUnit) < 75.GB) {
        log.warn "--read_aligner bwamem3 builds the ${params.host_db} decontamination index with ~75 GB RAM, but --max_memory is ${params.max_memory}: the build will likely fail. Raise --max_memory, pass a pre-built bwa-mem3 index with --custom_decontamination_index, or use --read_aligner bowtie2 (~6 GB to build)"
    }

    // A pre-built --custom_decontamination_index must match --read_aligner: a
    // Bowtie2 index cannot serve bwa-mem3 (a bwa-mem2 index can), nor the
    // reverse. Checked here, before any channel exists: raised later, while
    // INPUT_CHECK's samplesheet channel was starting, the error could be
    // reported as "Cannot find samplesheet file" instead
    if (flagOn(params.quality_control) && params.custom_decontamination_index) {
        if ( params.read_aligner == 'bwamem3' ) {
            def index_dir = file(params.custom_decontamination_index, checkIfExists: true)
            def index_files = index_dir.isDirectory() ? index_dir.listFiles().collect { f -> f.name } : []
            if ( !index_files.any { name -> name.endsWith('.bwt.2bit.64') } ) {
                def bowtie2_hint = index_files.any { name -> name ==~ /.*\.bt2l?$/ } ?
                    " It holds a Bowtie2 index: run with --read_aligner bowtie2 to keep using it, rebuild it with bwa-mem3 (or bwa-mem2), or drop --custom_decontamination_index and pass --custom_host_fasta." : ""
                error "--read_aligner bwamem3 needs --custom_decontamination_index to be a directory with a bwa-mem3/bwa-mem2 index (*.bwt.2bit.64): ${params.custom_decontamination_index}.${bowtie2_hint}"
            }
        } else {
            // The reverse mix-up: the index the pipeline publishes by
            // default (bwamem3_index/) passed under the bowtie2 fallback
            def index_dir = file(params.custom_decontamination_index, checkIfExists: true)
            def index_files = index_dir.isDirectory() ? index_dir.listFiles().collect { f -> f.name } : []
            if ( !index_files.any { name -> name ==~ /.*\.[0-9]\.bt2l?$/ } ) {
                def bwamem3_hint = index_files.any { name -> name.endsWith('.bwt.2bit.64') } ?
                    " It holds a bwa-mem3/bwa-mem2 index: drop --read_aligner bowtie2 to use it with the default aligner, or point the option at a Bowtie2 index." : ""
                error "--read_aligner bowtie2 needs --custom_decontamination_index to be a directory with a Bowtie2 index (*.bt2): ${params.custom_decontamination_index}.${bwamem3_hint}"
            }
        }
    }

    // Print run configuration
    log.info "  Run configuration:"
    log.info "  -------------------"
    log.info "  Input samplesheet    : ${params.input}"
    log.info "  Output directory     : ${params.output}"
    log.info "  Quality control      : ${params.quality_control}"
    log.info "  Read aligner         : ${params.read_aligner}"
    log.info "  Assembly mode        : ${params.assembly_mode}"
    log.info "  Taxonomic profiler   : ${params.taxonomic_profiler}${params.taxonomic_profiler == 'sylph' ? " (${params.custom_sylph_db ? 'custom database: ' : ''}${sylphDbLabel()})" : ''}"
    log.info "  Include binning      : ${params.include_binning}"
    log.info "  Binners              : ${binners_list.join(', ')}${binners_list.size() >= 2 ? ' (+ MetaWRAP refinement)' : ''}"
    log.info "  Read ARG prediction  : ${params.read_arg_prediction}"
    log.info "  RGI AMR prediction   : ${params.rgi_prediction}"
    log.info "  K-MARVEL prediction  : ${params.kmarvel_prediction}"
    log.info "  Contig tax and ARG   : ${params.contig_tax_and_arg}"
    log.info "  Contig functional    : ${params.contig_level_functional}"
    log.info "  MAG functional (Bakta): ${params.mag_level_functional}"
    log.info "  Read functional      : ${params.read_level_functional}${params.read_level_functional == 'superfocus' ? " (${params.superfocus_aligner})" : params.read_level_functional == 'humann' ? " (${params.custom_humann_db ? 'custom database' : params.humann_db})" : ''}"
    log.info "  MicrobeCensus        : ${params.microbecensus}"
    log.info "  dbCAN CAZy           : ${params.functional_cazy}"
    log.info ""

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
    // SUBWORKFLOW: Quality control and decontamination (single-pass phiX + host removal)
    //
    QC(
        ch_reads,
        PREPARE_DATABASES.out.decontamination_index
    )
    
    // In DSL2, process outputs can be referenced multiple times
    // No need to split channels - just use QC.out.reads directly in each consumer
    ch_clean_reads           = QC.out.reads
    ch_reads_report          = QC.out.report

    // Software provenance: every stage mixes its versions.yml files in here;
    // they are aggregated into pipeline_info/software_versions.yml at the end
    ch_versions = channel.empty()
        .mix(PREPARE_DATABASES.out.versions)
        .mix(QC.out.versions)

    //
    // SUBWORKFLOW: Taxonomic profiling
    //
    // Inputs of the run-level summary; each stays empty when its stage is off
    ch_summary_tax_report = channel.empty()
    ch_summary_tax_tables = channel.empty()
    ch_summary_sylph      = channel.empty()
    ch_summary_kmarvel    = channel.empty()
    ch_summary_kmarvel_db = channel.empty()
    ch_summary_rgi        = channel.empty()
    ch_summary_rgi_kmer   = channel.empty()

    if ( params.taxonomic_profiler != "none" ) {
        TAXONOMY(
            ch_clean_reads,
            ch_reads_report,
            PREPARE_DATABASES.out.kraken_db,
            PREPARE_DATABASES.out.sourmash_db,
            PREPARE_DATABASES.out.sylph_db,
            PREPARE_DATABASES.out.sylph_taxonomy
        )
        ch_versions = ch_versions.mix(TAXONOMY.out.versions)
        ch_summary_tax_report = TAXONOMY.out.reads_report
        ch_summary_tax_tables = TAXONOMY.out.tables
        ch_summary_sylph      = TAXONOMY.out.sylph_profiles
    }

    //
    // ARG PREDICTION IN READS
    //
    if ( flagOn(params.read_arg_prediction) ) {
        ARGS_OAP(ch_clean_reads)
        ch_args_oap = ARGS_OAP.out.args_oap_s1
        ch_argv_prediction = KARGVA(ch_clean_reads.combine(PREPARE_DATABASES.out.kargva_db))
        KARGA(ch_argv_prediction.kargva_reads.combine(PREPARE_DATABASES.out.karga_db))
        ch_arg_prediction = KARGA.out.kargva
        ARG_NORM_REPORT(
            ch_arg_prediction
                .concat(ch_argv_prediction.kargva_reports)
                .concat(ch_args_oap)
                .collect(sort: true)
        )
        ch_versions = ch_versions.mix(
            ARGS_OAP.out.versions.first(),
            KARGVA.out.versions.first(),
            KARGA.out.versions.first(),
            ARG_NORM_REPORT.out.versions
        )
    }

    //
    // RGI AMR PREDICTION IN READS
    //
    if ( flagOn(params.rgi_prediction) ) {
        // Store rgi_card_db to allow reuse
        ch_rgi_db = PREPARE_DATABASES.out.rgi_card_db
            .ifEmpty { error "ERROR: RGI database is empty. Ensure params.rgi_prediction is enabled and a valid CARD database is configured." }
        
        // RGI bwt: Align reads to CARD AMR alleles
        ch_rgi_bwt = RGI_BWT(
            ch_clean_reads,
            ch_rgi_db.collect()
        )
        ch_summary_rgi = ch_rgi_bwt.allele_mapping.map { _meta, file -> file }
        
        // RGI kmer: Pathogen-of-origin prediction. The k-mer database comes
        // with WildCARD, so the downloaded CARD-only database has none; a
        // --custom_rgi_card_db must have it loaded (RGI_KMER checks)
        ch_rgi_kmer_json = channel.value([])
        if ( params.custom_rgi_card_db || flagOn(params.rgi_include_wildcard) ) {
            ch_rgi_kmer = RGI_KMER(
                ch_rgi_bwt.bam,
                ch_rgi_db.collect()
            )
            ch_rgi_kmer_json = ch_rgi_kmer.kmer_json.map { _meta, file -> file }.collect(sort: true).ifEmpty([])
            ch_summary_rgi_kmer = ch_rgi_kmer.gene_predictions.map { _meta, file -> file }
            ch_versions = ch_versions.mix(ch_rgi_kmer.versions.first())
        } else {
            log.warn "RGI pathogen of origin (rgi kmer_query) is skipped: its k-mer database comes with WildCARD (--rgi_include_wildcard true)"
        }

        // Generate summary report
        RGI_REPORT(
            ch_rgi_bwt.allele_mapping.map { _meta, file -> file }.collect(sort: true),
            ch_rgi_bwt.gene_mapping.map { _meta, file -> file }.collect(sort: true),
            ch_rgi_kmer_json
        )
        ch_versions = ch_versions.mix(
            ch_rgi_bwt.versions.first(),
            RGI_REPORT.out.versions
        )
    }

    //
    // K-MARVEL AMR PREDICTION IN READS
    //
    if ( flagOn(params.kmarvel_prediction) ) {
        ch_kmarvel_db = PREPARE_DATABASES.out.kmarvel_db
            .ifEmpty { error "ERROR: K-MARVEL database is empty. Ensure params.kmarvel_prediction is enabled and a valid CARD database is configured." }

        // Assemble CARD ARGs from each sample's reads
        KMARVEL(
            ch_clean_reads,
            ch_kmarvel_db.collect()
        )

        // Merge the per-sample tables into study-level matrices
        KMARVEL_REPORT(
            KMARVEL.out.tsv.map { _meta, file -> file }.collect(sort: true)
        )
        ch_versions = ch_versions.mix(
            KMARVEL.out.versions.first(),
            KMARVEL_REPORT.out.versions
        )
        ch_summary_kmarvel    = KMARVEL.out.tsv.map { _meta, file -> file }
        ch_summary_kmarvel_db = ch_kmarvel_db
    }

    //
    // SUBWORKFLOW: Assembly
    //
    ch_contigs_meta  = channel.empty()
    ch_bam_meta      = channel.empty()
    ch_counting_bam  = channel.empty()
    ch_refined_bins  = channel.empty()

    if ( params.assembly_mode != "none" ) {
        ASSEMBLY(
            ch_clean_reads
        )

        ch_contigs_meta = ASSEMBLY.out.contigs_meta
        ch_bam_meta     = ASSEMBLY.out.bam_meta
        ch_counting_bam = ASSEMBLY.out.counting_bam
        ch_versions     = ch_versions.mix(ASSEMBLY.out.versions)

        //
        // SUBWORKFLOW: Binning
        //
        if ( flagOn(params.include_binning) ) {
            BINNING(
                ASSEMBLY.out.bam,
                PREPARE_DATABASES.out.checkm2_db
                    .ifEmpty { error "ERROR: CheckM2 database is empty. Ensure params.include_binning is enabled and a valid CheckM2 database is configured." },
                PREPARE_DATABASES.out.gtdbtk_db
                    .ifEmpty { error "ERROR: GTDB-Tk database is empty. Ensure params.include_binning is enabled and a valid GTDB-Tk database is configured." },
                ch_clean_reads
            )
            ch_refined_bins = BINNING.out.refined_bins
            ch_versions = ch_versions.mix(BINNING.out.versions)
        }
    }

    //
    // SHARED GENE CALLING ON CONTIGS
    // One Pyrodigal pass (metagenome mode, both assembly modes) whose FAA feeds
    // DeepARG and whose FAA/GFF feed the functional annotation branch (Q9).
    //
    ch_contig_proteins  = channel.empty()
    ch_contig_genes_gff = channel.empty()

    if ( (flagOn(params.contig_tax_and_arg) || flagOn(params.contig_level_functional)) && params.assembly_mode != "none" ) {
        //
        // Run nf-core PYRODIGAL for contig ORF prediction
        // nf-core PYRODIGAL signature:
        //   input:  tuple val(meta), path(fasta) + val(output_format)
        //   output: tuple val(meta), path("*.faa.gz"), emit: faa
        //           tuple val(meta), path("*.gff.gz"), emit: annotations
        //
        PYRODIGAL(
            ch_contigs_meta,
            "gff"  // output_format
        )
        ch_contig_proteins  = PYRODIGAL.out.faa
        ch_contig_genes_gff = PYRODIGAL.out.annotations
        ch_versions = ch_versions.mix(PYRODIGAL.out.versions.first())
    }

    //
    // MODULE: MicrobeCensus average genome size on host-removed reads
    // (design doc Section 4.7). Sits outside the branch subworkflows because
    // both the contig and the read branch (T8) consume its output; the read
    // branch needs no assembly.
    // Failure is non-fatal (errorStrategy in config/modules.config): a failed
    // sample emits nothing here and falls back to TPM-only in aggregation.
    def run_microbecensus = flagOn(params.microbecensus)
    ch_ags = channel.empty()
    def run_contig_functional = flagOn(params.contig_level_functional) && params.assembly_mode != "none"
    def run_read_functional   = params.read_level_functional != "none"
    if ( run_microbecensus && (run_contig_functional || run_read_functional) ) {
        MICROBECENSUS(ch_clean_reads)
        ch_ags = MICROBECENSUS.out.ags
        ch_versions = ch_versions.mix(MICROBECENSUS.out.versions.first())
    }

    //
    // SUBWORKFLOW: Functional annotation — contig branch (eggNOG-mapper +
    // run_dbcan CAZy + featureCounts gene quantification) and/or MAG branch
    // (Bakta on refined bins). The branches are independent; each DB channel
    // is only demanded (.ifEmpty error) when its branch is on
    //
    if ( (flagOn(params.contig_level_functional) || flagOn(params.mag_level_functional)) && params.assembly_mode != "none" ) {
        // With functional_cazy off the dbCAN DB channel is legitimately empty
        // and the subworkflow never consumes it
        def run_functional_cazy = flagOn(params.contig_level_functional) && flagOn(params.functional_cazy)
        ch_dbcan_db = run_functional_cazy
            ? PREPARE_DATABASES.out.dbcan_db
                .ifEmpty { error "ERROR: dbCAN database is empty. Ensure params.functional_cazy is enabled and a valid dbCAN database is configured." }
            : channel.empty()
        ch_eggnog_db = flagOn(params.contig_level_functional)
            ? PREPARE_DATABASES.out.eggnog_db
                .ifEmpty { error "ERROR: eggNOG database is empty. Ensure params.contig_level_functional is enabled and a valid eggNOG database is configured." }
            : channel.empty()
        ch_cog_def = flagOn(params.contig_level_functional)
            ? PREPARE_DATABASES.out.cog_def
                .ifEmpty { error "ERROR: COG definitions table is empty. --contig_level_functional (and --read_level_functional woltka) need a valid COG table (--cog_db, or --custom_cog_db)." }
            : channel.empty()
        ch_bakta_db = flagOn(params.mag_level_functional)
            ? PREPARE_DATABASES.out.bakta_db
                .ifEmpty { error "ERROR: Bakta database is empty. Ensure params.mag_level_functional is enabled and a valid Bakta database is configured." }
            : channel.empty()
        FUNCTIONAL_ANNOTATION(
            ch_contig_proteins,
            ch_eggnog_db,
            ch_dbcan_db,
            ch_contig_genes_gff,
            ch_counting_bam,
            ch_ags,
            ch_refined_bins,
            ch_bakta_db,
            ch_cog_def
        )
        ch_versions = ch_versions.mix(FUNCTIONAL_ANNOTATION.out.versions)
    }

    //
    // SUBWORKFLOW: Read-level functional profiling (design doc Section 4.6),
    // one backend per run, independent of assembly; its read_* tables are
    // reported separately from the contig branch's
    //
    if ( run_read_functional ) {
        // the selected backend's database (only that one is provisioned)
        def ch_read_db = params.read_level_functional == 'woltka' ?
            PREPARE_DATABASES.out.woltka_db
                .ifEmpty { error "ERROR: Woltka (WoLr2) database is empty. Ensure --read_level_functional woltka is set and a valid WoLr2 database is configured (--custom_woltka_db)." } :
            params.read_level_functional == 'humann' ?
            PREPARE_DATABASES.out.humann_db
                .ifEmpty { error "ERROR: HUMAnN database is empty. Ensure --read_level_functional humann is set and a valid HUMAnN database root is configured (--custom_humann_db: chocophlan/, uniref/, utility_mapping/, metaphlan/)." } :
            PREPARE_DATABASES.out.superfocus_db
                .ifEmpty { error "ERROR: SUPER-FOCUS database is empty. Ensure --read_level_functional superfocus is set and a valid SUPER-FOCUS database root is configured (--custom_superfocus_db, the directory containing db/)." }
        // the woltka backend maps WoLr2 COG ids to categories with the
        // contig branch's NCBI COG table (design doc Q14)
        def ch_read_cog_def = params.read_level_functional == 'woltka' ?
            PREPARE_DATABASES.out.cog_def
                .ifEmpty { error "ERROR: COG definitions table is empty. --read_level_functional woltka needs a valid COG table (--cog_db, or --custom_cog_db)." } :
            channel.empty()
        READ_FUNCTIONAL(
            ch_clean_reads,
            ch_read_db,
            ch_ags,
            ch_read_cog_def
        )
        ch_versions = ch_versions.mix(READ_FUNCTIONAL.out.versions)
    }

    //
    // CONTIG-LEVEL TAXONOMY AND ARG PREDICTION
    //
    if ( flagOn(params.contig_tax_and_arg) && params.assembly_mode != "none" ) {
        // The list wrap keeps a multi-file BLAST DB as ONE tuple element
        // (path(nt_db)) instead of flattening it into the tuple, which staged
        // only the first DB file into NT_BLASTN
        NT_BLASTN(ch_contigs_meta.combine(PREPARE_DATABASES.out.blast_db
            .ifEmpty { error "ERROR: BLAST database is empty. Ensure params.contig_tax_and_arg is enabled and a valid BLAST database is configured." }
            .map { db_files -> [db_files] }))
        ch_nt_blastn = NT_BLASTN.out.megablast_to_blob
        //
        // Run nf-core SAMTOOLS_INDEX for BAM indexing
        // nf-core SAMTOOLS_INDEX signature:
        //   input:  tuple val(meta), path(input)
        //   output: tuple val(meta), path("*.bai"), emit: bai
        //
        NFCORE_SAMTOOLS_INDEX(ch_bam_meta)
        
        // Join BAM with its index for BLOBTOOLS compatibility
        // BLOBTOOLS expects: tuple val(meta), path(bam), path(bam_bai)
        ch_index_bam = ch_bam_meta
            .join(NFCORE_SAMTOOLS_INDEX.out.bai)
        // Wrap the collected taxdump files in a list so they arrive as ONE
        // tuple element (path(tax_files)) instead of being flattened into the
        // tuple — flattened, only the first file was staged into BLOBTOOLS
        ch_blob_table = BLOBTOOLS(
            ch_nt_blastn
                .join(ch_index_bam)
                .combine(PREPARE_DATABASES.out.taxdump.collect().map { files -> [files] })
        )
        BLOBPLOT(ch_blob_table.only_blob.collect(sort: true))

        // Contig proteins come from the shared PYRODIGAL step above
        ch_contig_args = DEEPARG_CONTIGS(ch_contig_proteins.combine(PREPARE_DATABASES.out.deeparg_db
            .ifEmpty { error "ERROR: DeepARG database is empty. Ensure params.contig_tax_and_arg is enabled and a valid DeepARG database is configured." }))
        ch_arg_contig_data = ARG_CONTIG_LEVEL_REPORT(
            ch_contig_args.only_deeparg
                .concat(ch_blob_table.only_blob)
                .collect(sort: true)
        )
        ARG_BLOBPLOT(ch_arg_contig_data.arg_reports)

        ch_versions = ch_versions.mix(
            NT_BLASTN.out.versions.first(),
            NFCORE_SAMTOOLS_INDEX.out.versions.first(),
            BLOBTOOLS.out.versions.first(),
            BLOBPLOT.out.versions,
            DEEPARG_CONTIGS.out.versions.first(),
            ARG_CONTIG_LEVEL_REPORT.out.versions,
            ARG_BLOBPLOT.out.versions
        )
    }

    //
    // ARG PREDICTION IN BINS AND CLUSTERING
    //
    if ( flagOn(params.arg_bin_clustering) && flagOn(params.include_binning) ) {
        PRODIGAL_BINS(ch_refined_bins)
        ch_raw_orfs = PRODIGAL_BINS.out.prodigal_bins
        DEEPARG_BINS(ch_raw_orfs.combine(PREPARE_DATABASES.out.deeparg_db
            .ifEmpty { error "ERROR: DeepARG database is empty. Ensure params.arg_bin_clustering is enabled and a valid DeepARG database is configured." }))
        ch_deeparg = DEEPARG_BINS.out.deeparg_bins
        ARG_FASTA_FORMATTER(ch_raw_orfs.join(ch_deeparg))
        ch_arg_fasta = ARG_FASTA_FORMATTER.out.arg_reports
        CLUSTERING(ch_arg_fasta.collect(sort: true))

        ch_versions = ch_versions.mix(
            PRODIGAL_BINS.out.versions.first(),
            DEEPARG_BINS.out.versions.first(),
            ARG_FASTA_FORMATTER.out.versions.first(),
            CLUSTERING.out.versions
        )
    }

    //
    // Aggregate software versions -> pipeline_info/software_versions.yml
    // Per-sample tasks of one process write byte-identical versions.yml files,
    // so content-level unique() collapses them to one block per process.
    // stripIndent() normalizes the varying heredoc indentation across modules.
    //
    ch_software_versions = ch_versions
        .map { yml -> yml.text.stripIndent() }
        .mix( channel.of(
            ( "\"${workflow.manifest.name}\":\n" +
              "    pipeline: ${workflow.manifest.version}\n" +
              "    nextflow: ${nextflow.version}\n" ).toString()
        ) )
        .unique()
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.output}/pipeline_info",
            sort: true,
            newLine: false
        )

    //
    // MODULE: Run-level summary (00_summary/): QC, taxonomy and read-level
    // ARG (K-MARVEL, RGI) tables, Excel workbook, figures and HTML report.
    // It reads the aggregated software_versions.yml, so its own versions.yml
    // is not part of that file.
    //
    if ( flagOn(params.results_summary) ) {
        // the profiler's Reads_report.csv adds the classified fraction
        // (TAXONOMY_REPORT's output is already a value channel)
        ch_summary_reads_report = params.taxonomic_profiler != 'none'
            ? ch_summary_tax_report
            : ch_reads_report.flatten().filter { report -> report.name.endsWith('.csv') }.first()
        RESULTS_SUMMARY(
            ch_summary_reads_report,
            QC.out.fastp_json.toSortedList { a, b -> a.name <=> b.name },
            ch_summary_tax_tables.toSortedList { a, b -> a.name <=> b.name },
            ch_summary_sylph.toSortedList { a, b -> a.name <=> b.name },
            ch_summary_kmarvel.toSortedList { a, b -> a.name <=> b.name },
            ch_summary_kmarvel_db.toList(),
            ch_summary_rgi.toSortedList { a, b -> a.name <=> b.name },
            ch_summary_rgi_kmer.toSortedList { a, b -> a.name <=> b.name },
            ch_software_versions,
            params.sample_metadata ? file(params.sample_metadata, checkIfExists: true) : [],
            params.card_prevalence ? file(params.card_prevalence, checkIfExists: true) : [],
            file(params.arg_priority_table, checkIfExists: true),
            file(params.summary_watchlist, checkIfExists: true)
        )
    }

    //
    // COMPLETION HANDLERS
    // Registered at the end of the workflow body (strict syntax disallows
    // top-level handlers): a validation error() above throws before these run,
    // so they do not fire on startup failures — same behavior as before.
    //
    workflow.onComplete = {
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

    workflow.onError = {
        log.error "Pipeline failed. Check error message above or in ${params.output}/pipeline_info/"
    }
}
