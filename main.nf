#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gene2dis/BugBuster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/gene2dis/BugBuster
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE LOGO AND INFO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

// Input validation
include { INPUT_CHECK } from './subworkflows/local/input_check'

	// FORMATEO Y/O DESCARGA DE BASES DE DATOS

include {FORMAT_SM_DB} from "./modules/format_db/main"
include {FORMAT_KRAKEN_DB} from "./modules/format_db/main"
include {FORMAT_BOWTIE_INDEX as FORMAT_BOWTIE_HOST} from "./modules/format_db/main"
include {FORMAT_NT_BLAST_DB} from "./modules/format_db/main"
include {FORMAT_TAXDUMP_FILES} from "./modules/format_db/main"
include {DOWNLOAD_DEEPARG_DB} from "./modules/format_db/main"
include {FORMAT_CHECKM2_DB} from "./modules/format_db/main"
include {BUILD_PHIX_BOWTIE2_INDEX} from "./modules/format_db/main"
include {DOWNLOAD_GTDBTK_DB} from "./modules/format_db/main"

	// SEQUENCE FILTERING

include {FASTP} from "./modules/fastp/main"
include {QFILTER} from "./modules/qfilter/main"
include {COUNT_READS} from "./modules/count_reads/main"

include {BOWTIE2 as BOWTIE2_HOST} from "./modules/bowtie2/main"
include {BOWTIE2 as BOWTIE2_PHYX} from "./modules/bowtie2/main"

	// TAXONOMIC PREDICTION IN READS

include {KRAKEN2 as KRAKEN2} from "./modules/kraken2/main"
include {BRACKEN as BRACKEN} from "./modules/bracken/main"
include {KRAKEN_BIOM as KRAKEN_BIOM} from "./modules/kraken_biom/main"
include {KRAKEN_TO_PHYLOSEQ} from "./modules/kraken_to_phyloseq/main"


include {SOURMASH} from "./modules/sourmash/main"
include {SOURMASH_TO_PHYLOSEQ} from "./modules/sourmash_to_phyloseq/main"

	// ASSEMBLY AND CO-ASSEMBLY

include {MEGAHIT} from "./modules/megahit/main"
include {MEGAHIT_COASSEMBLY} from "./modules/megahit/main"

include {BBMAP} from "./modules/bbmap/main"
include {BBMAP_COASSEMBLY} from "./modules/bbmap/main"

	// BINNING 

include {BOWTIE2_SAMTOOLS} from "./modules/bowtie2_samtools/main"
include {BOWTIE2_SAMTOOLS_COASSEMBLY} from "./modules/bowtie2_samtools/main"
include {BOWTIE2_SAMTOOLS_DEPTH} from "./modules/bowtie2_samtools/main"

include {CALCULATE_DEPTH} from "./modules/calculate_depth/main"
include {CALCULATE_DEPTH_COASSEMBLY} from "./modules/calculate_depth/main"

include {SEMIBIN} from "./modules/semibin/main"
include {SEMIBIN_COASSEMBLY} from "./modules/semibin/main"

include {METABAT2} from "./modules/metabat2/main"
include {METABAT2_COASSEMBLY} from "./modules/metabat2/main"

include {AUTOMETA} from "./modules/autometa/main"
include {AUTOMETA_COASSEMBLY} from "./modules/autometa/main"

include {COMEBIN} from "./modules/comebin/main"
include {COMEBIN_COASSEMBLY} from "./modules/comebin/main"

include {METAWRAP} from "./modules/metawrap/main"
include {METAWRAP_COASSEMBLY} from "./modules/metawrap/main"

	// BIN QUALITY PREDICTION

include {CHECKM2} from "./modules/checkm2/main"
include {CHECKM2_COASSEMBLY} from "./modules/checkm2/main"

include {BIN_QUALITY_REPORT} from "./modules/bin_quality_report/main"
include {BIN_TAX_REPORT} from "./modules/bin_tax_report/main"
include {BIN_SUMMARY} from "./modules/bin_summary/main"

	// FUNCTIONAL ANNOTATION

include {METACERBERUS_CONTIGS} from "./modules/metacerberus/main"

	// TAXONOMIC PREDICTION IN CONTIGS

include {NT_BLASTN} from "./modules/nt_blastn/main"
include {BLOBTOOLS} from "./modules/blobtools/main"
include {SAMTOOLS_INDEX} from "./modules/samtools_index/main"
include {BLOBPLOT} from "./modules/blobplot/main"

	// READ REPORTS

include {READS_REPORT} from "./modules/reads_report/main"
include {TAX_REPORT_KRAKEN2} from "./modules/tax_report_kraken2/main"
include {TAX_REPORT_SOURMASH} from "./modules/tax_report_sourmash/main"

	// BIN TAXONOMIC PREDICTION

include {GTDB_TK} from "./modules/gtdb-tk/main"
include {GTDB_TK_COASSEMBLY} from "./modules/gtdb-tk/main"

	// COVERAGE CALCULATION PER BIN

include {BEDTOOLS} from "./modules/bedtools/main"

	// ORF PREDICTION IN CONTIGS AND BINS

include {PRODIGAL_BINS} from "./modules/prodigal/main"
include {PRODIGAL_CONTIGS} from "./modules/prodigal/main"

	// ARG PREDICTION IN READS

include {KARGVA} from "./modules/kargva/main"
include {KARGA} from "./modules/karga/main"
include {ARGS_OAP} from "./modules/args_oap/main"
include {ARG_NORM_REPORT} from "./modules/arg_norm_report/main"

	// ARG PREDICTION IN CONTIGS AND BINS

include {DEEPARG_BINS} from "./modules/deeparg/main"
include {DEEPARG_CONTIGS} from "./modules/deeparg/main"
include {ARG_CONTIG_LEVEL_REPORT} from "./modules/arg_contig_level_report/main"
include {ARG_FASTA_FORMATTER} from "./modules/arg_fasta_formatter/main"
include {CLUSTERING} from "./modules/clustering/main"
include {ARG_BLOBPLOT} from "./modules/arg_blobplot/main"

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
    // DATABASE SETUP: Taxonomic profiler databases
    //
    switch ( params.taxonomic_profiler ) {

       case "kraken2":

       if ( params.custom_kraken_db ) {
             ch_kraken_ref_db_formated = Channel.fromPath(params.custom_kraken_db)
       } else {
             ch_kraken_ref_db = Channel.fromList(params.kraken_ref_db[params.kraken2_db]["file"])
             ch_kraken_ref_db_formated = FORMAT_KRAKEN_DB(ch_kraken_ref_db)
       }

       break;

       case "sourmash":

       if ( params.custom_sourmash_db ) {
             ch_sourmash_db = Channel.fromList(params.custom_sourmash_db).map { file(it) }.collect()
       } else {
             ch_sourmash_db = Channel.fromList(params.sourmash_ref_db[params.sourmash_db]["file"]).map { file(it) }
       }

       break;

       case "none":

       break;

       default:
       exit 1, "ERROR: Not valid taxonomic profiler, please check nextflow.config file"
    }

    if ( params.read_arg_prediction ) {
       if ( params.custom_karga_db ) {
             ch_karga_db = Channel.from(file(params.custom_karga_db))
       } else {
             ch_karga_db = Channel.fromList(params.karga_ref_db[params.karga_db]["file"]).map { file(it) }
       }
    

    if ( params.custom_kargva_db ) {
             ch_kargva_db = Channel.from(file(params.custom_kargva_db))
       } else {
             ch_kargva_db = Channel.fromList(params.kargva_ref_db[params.kargva_db]["file"]).map { file(it) }
       }
    }

    if ( params.quality_control ) {

       if ( params.custom_phiX_index ) {
            ch_bowtie_phiX_index_formated = Channel.fromPath(params.custom_phiX_index)
       } else {
            ch_bowtie_phiX_index = Channel.fromList(params.bowtie_ref_genomes_for_build[params.phiX_index]["file"]).map { file(it) }
            ch_bowtie_phiX_index_formated = BUILD_PHIX_BOWTIE2_INDEX(ch_bowtie_phiX_index)
       }

       if ( params.custom_bowtie_host_index ) {
            ch_bowtie_host_index_formated = Channel.fromPath(params.custom_bowtie_host_index)
       } else {
            ch_bowtie_host_index = Channel.fromList(params.bowtie_ref_host_index[params.host_db]["file"])
            ch_bowtie_host_index_formated = FORMAT_BOWTIE_HOST(ch_bowtie_host_index)
       }
    }

    if ( params.contig_tax_and_arg ) {

       if ( params.custom_deeparg_db ) {
           ch_deeparg_db = Channel.fromPath(params.custom_deeparg_db)
       } else {
           ch_deeparg_db = DOWNLOAD_DEEPARG_DB()
       }

       if ( params.custom_blast_db ) {
           ch_blast_db_formated = Channel.fromPath(params.custom_blast_db)
       } else {
           ch_blast_db =  Channel.fromList(params.blast_ref_db[params.blast_db]["url"])
           ch_blast_db_formated = FORMAT_NT_BLAST_DB(ch_blast_db)
       }

       if ( params.custom_taxdump_files ) {
           ch_taxdump_files_formated = Channel.from(path(params.custom_taxdump_files))
       } else {
           ch_taxdump_files = Channel.fromList(params.taxonomy_files[params.taxdump_files]["url"])
           ch_taxdump_files_formated = FORMAT_TAXDUMP_FILES(ch_taxdump_files)
       }
    }
    
    if ( params.include_binning ) {

       if ( params.custom_gtdbtk_db ) {
           ch_gtdbtk_db_formated = Channel.fromPath(params.custom_gtdbtk_db)
       } else {
           ch_gtdbtk_db = Channel.fromList(params.gtdbtk_ref_db[params.gtdbtk_db]["url"])
           ch_gtdbtk_db_formated = DOWNLOAD_GTDBTK_DB(ch_gtdbtk_db)
       }

       if ( params.custom_checkm2_db ) {
           ch_checkm2_db_formated = Channel.fromPath(params.custom_checkm2_db)
       } else {
           ch_checkm2_db = Channel.fromList(params.checkm2_ref_db[params.checkm2_db]["url"])
           ch_checkm2_db_formated = FORMAT_CHECKM2_DB(ch_checkm2_db)
       }

    }

    //
    // QUALITY CONTROL: Read filtering and host decontamination
    //
    if ( params.quality_control ) {

        // Read quality filtering with Fastp
        ch_fastp_reads = FASTP(ch_reads)

        // Extract and format QC reports
        ch_fastp_reads_report = QFILTER(ch_fastp_reads.fastq)

        // Filter samples by minimum read count
        ch_fastp_reads_num = ch_fastp_reads_report.qfilter.map { tuple ->
            def after_reads = tuple[2].text
            if ( Integer.parseInt(after_reads) >= params.min_read_sample ) {
                def meta = tuple[0]
                def reads = tuple[1]
                return [meta, reads]
            } else {
                return null
            }
        }.filter { it != null }

        // Remove PhiX and host contamination
        ch_phyx_clean_reads = BOWTIE2_PHYX(ch_fastp_reads_num.combine(ch_bowtie_phiX_index_formated), "phiX")
        ch_host_clean_reads = BOWTIE2_HOST(ch_phyx_clean_reads.reads.combine(ch_bowtie_host_index_formated), "host")

        // Collect read reports
        ch_reads_report = READS_REPORT(
            ch_phyx_clean_reads.report
                .concat(ch_host_clean_reads.report)
                .concat(ch_fastp_reads_report.reads_report)
                .collect(),
            "phiX host"
        )
    } else {
        ch_host_clean_reads = COUNT_READS(ch_reads)

        ch_reads_report = READS_REPORT(
            ch_host_clean_reads.reads_report.collect(), 
            "none"
        )
    }

    switch ( params.taxonomic_profiler ) {
       case "kraken2":

            ch_k2_taxonomy = KRAKEN2(ch_host_clean_reads.reads.combine(ch_kraken_ref_db_formated), params.kraken_db_used)

            TAX_REPORT_KRAKEN2(ch_k2_taxonomy.report
                           .concat(ch_reads_report)
                           .collect())

        // ABUNDANCE ESTIMATION - requires specifying the database path,
        // an output alias, minimum read size and taxonomic level.
        // If this step fails, try changing the taxonomic level.

            ch_taxonomy_estimation = BRACKEN(ch_k2_taxonomy.kraken.combine(ch_kraken_ref_db_formated), params.kraken_db_used)

        // GENERATE .BIOM FILE

            ch_kraken_biom = KRAKEN_BIOM(ch_taxonomy_estimation.collect(), params.kraken_db_used)

        // GENERATE PHYLOSEQ OBJECT AND ABUNDANCE PLOTS

            KRAKEN_TO_PHYLOSEQ(ch_kraken_biom)

       break;

       case "sourmash":

            ch_sm_taxonomy = SOURMASH(ch_host_clean_reads.reads.combine(ch_sourmash_db.collect()), params.sourmash_db_name, params.sourmash_tax_rank)

            TAX_REPORT_SOURMASH(ch_sm_taxonomy.report
                           .concat(ch_reads_report)
                           .collect())

            SOURMASH_TO_PHYLOSEQ(ch_sm_taxonomy.sourmash_gather.collect(), params.sourmash_db_name)


       break;

       case "none":

       break;

       default:

       exit 1, "ERROR: Not valid taxonomic profiler, please check nextflow.config file"
    }

	// ARG PREDICTION IN READS

    if ( params.read_arg_prediction) {

          ch_args_oap = ARGS_OAP(ch_host_clean_reads.reads)
          ch_argv_prediction = KARGVA(ch_host_clean_reads.reads.combine(ch_kargva_db))
          ch_arg_prediction = KARGA(ch_argv_prediction.kargva_reads.combine(ch_karga_db))
          ARG_NORM_REPORT(ch_arg_prediction.concat(ch_argv_prediction.kargva_reports).concat(ch_args_oap).collect())

    }

    switch ( params.assembly_mode ) {
       case "coassembly": 

        // READ ASSEMBLY, CONTIG FILTERING AND DEPTH CALCULATION

         ch_assembly_co = MEGAHIT_COASSEMBLY(ch_host_clean_reads.reads_coassembly.collect())
         ch_filtered_contigs_co = BBMAP_COASSEMBLY(ch_assembly_co)

         if ( params.include_binning || params.contig_tax_and_arg ) { 

              ch_bowtie2_samtools = BOWTIE2_SAMTOOLS_COASSEMBLY(ch_host_clean_reads.reads.combine(ch_filtered_contigs_co.contigs)) 
              ch_contigs_blastn = ch_filtered_contigs_co.contigs.map{
                                        def meta = [:]
                                            meta.id = "coassembly"
                                        def contigs = it
                                        return [meta, contigs]
              }

              ch_bam_index = ch_bowtie2_samtools.map{
                                        def meta = [:]
                                            meta.id = "coassembly"
                                        def bam = it
                                        return [meta, bam]
              }
         }

         if ( params.include_binning ) {

              ch_depth = CALCULATE_DEPTH_COASSEMBLY(ch_bowtie2_samtools.collect())

        // BINNING AND BIN REFINEMENT

              ch_metabat2 = METABAT2_COASSEMBLY(ch_filtered_contigs_co.contigs.combine(ch_depth))
              ch_semibin = SEMIBIN_COASSEMBLY(ch_bowtie2_samtools.collect().combine(ch_filtered_contigs_co.contigs))
              ch_comebin = COMEBIN_COASSEMBLY(ch_bowtie2_samtools.collect().combine(ch_filtered_contigs_co.contigs))

              ch_all_bins = ch_metabat2.combine(ch_semibin).combine(ch_comebin)
              ch_metawrap_co = METAWRAP_COASSEMBLY(ch_all_bins)

        // QUALITY ESTIMATION, TAXONOMIC ASSIGNMENT AND REPORT GENERATION

              ch_checkm = CHECKM2_COASSEMBLY(ch_all_bins.combine(ch_metawrap_co).combine(ch_checkm2_db_formated))
              ch_gtdb_tk = GTDB_TK_COASSEMBLY(ch_metawrap_co.combine(ch_gtdbtk_db_formated))
              ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(ch_host_clean_reads.reads.combine(ch_metawrap_co))
              ch_bin_cov = BEDTOOLS(ch_bin_depth)

              BIN_SUMMARY(ch_bin_cov.collect().combine(ch_gtdb_tk.report.collect()).combine(ch_checkm.metawrap_report.collect()))
         }

         break;
         
         case "assembly":

        // READ ASSEMBLY, CONTIG FILTERING AND DEPTH CALCULATION

         ch_assembly = MEGAHIT(ch_host_clean_reads.reads)
         ch_filtered_contigs = BBMAP(ch_assembly)

        // BINNING AND BIN REFINEMENT

         if ( params.include_binning || params.contig_tax_and_arg ) {
 
             ch_bowtie2_samtools = BOWTIE2_SAMTOOLS(ch_filtered_contigs.bbmap) 
             ch_contigs_blastn = ch_filtered_contigs.bbmap_contigs
             ch_bam_index = ch_bowtie2_samtools.only_bam

         }
                  
         if ( params.include_binning ) {

             ch_depth = CALCULATE_DEPTH(ch_bowtie2_samtools.contigs_and_bam)
             ch_metabat2 = METABAT2(ch_depth)
             ch_semibin = SEMIBIN(ch_bowtie2_samtools.contigs_and_bam)
             ch_comebin = COMEBIN(ch_bowtie2_samtools.contigs_and_bam)

             ch_all_bins = ch_metabat2.join(ch_semibin).join(ch_comebin)
             ch_metawrap = METAWRAP(ch_all_bins)

        // QUALITY ESTIMATION AND TAXONOMIC ASSIGNMENT

             ch_checkm = CHECKM2(ch_all_bins.join(ch_metawrap).combine(ch_checkm2_db_formated))
             ch_gtdb_tk = GTDB_TK(ch_metawrap.combine(ch_gtdbtk_db_formated))

             BIN_QUALITY_REPORT(ch_checkm.all_reports.collect())
             BIN_TAX_REPORT(ch_gtdb_tk.report.collect())
        }

         if ( params.contig_level_metacerberus) {
             METACERBERUS_CONTIGS(ch_filtered_contigs.bbmap_contigs)
         }  

        break;

        case "none":

        break;

        default:
         
        exit 1, "ERROR: Not valid assembly mode, please check nextflow.config file"
    }

    if ( params.contig_tax_and_arg ) {

             ch_nt_blastn = NT_BLASTN(ch_contigs_blastn.combine(ch_blast_db_formated))
             ch_index_bam = SAMTOOLS_INDEX(ch_bam_index)
             ch_blob_table = BLOBTOOLS(ch_nt_blastn.join(ch_index_bam).combine(ch_taxdump_files_formated.collect()))
             BLOBPLOT(ch_blob_table.only_blob.collect())

             ch_contig_proteins = PRODIGAL_CONTIGS(ch_contigs_blastn)
             ch_contig_args = DEEPARG_CONTIGS(ch_contig_proteins.combine(ch_deeparg_db))
             ch_arg_contig_data = ARG_CONTIG_LEVEL_REPORT(ch_contig_args.only_deeparg.concat(ch_blob_table.only_blob).collect())
             ARG_BLOBPLOT(ch_arg_contig_data)

    }

	// ARG PREDICTION IN BINS AND CLUSTERING

    if ( params.arg_bin_clustering) {

        ch_raw_orfs = PRODIGAL_BINS(ch_metawrap)
        ch_deeparg = DEEPARG_BINS(ch_raw_orfs.combine(ch_deeparg_db))     
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
