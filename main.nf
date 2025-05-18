// DSL-2 enabled
nextflow.enable.dsl = 2

/*
IMPORT MODULES
*/

logo='''
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
'''

println logo 

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

	// FILTRADO DE SECUENCIAS

include {FASTP} from "./modules/fastp/main"
include {QFILTER} from "./modules/qfilter/main"
include {COUNT_READS} from "./modules/count_reads/main"

include {BOWTIE2 as BOWTIE2_HOST} from "./modules/bowtie2/main"
include {BOWTIE2 as BOWTIE2_PHYX} from "./modules/bowtie2/main"

	// PREDICCION TAXONOMICA EN READS

include {KRAKEN2 as KRAKEN2} from "./modules/kraken2/main"
include {BRACKEN as BRACKEN} from "./modules/bracken/main"
include {KRAKEN_BIOM as KRAKEN_BIOM} from "./modules/kraken_biom/main"
include {KRAKEN_TO_PHYLOSEQ} from "./modules/kraken_to_phyloseq/main"


include {SOURMASH} from "./modules/sourmash/main"
include {SOURMASH_TO_PHYLOSEQ} from "./modules/sourmash_to_phyloseq/main"

	// ENSAMBLE Y COENSAMBLE 

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

	// PREDICCION DE CALIDAD EN BINS

include {CHECKM2} from "./modules/checkm2/main"
include {CHECKM2_COASSEMBLY} from "./modules/checkm2/main"

include {BIN_QUALITY_REPORT} from "./modules/bin_quality_report/main"
include {BIN_TAX_REPORT} from "./modules/bin_tax_report/main"
include {BIN_SUMMARY} from "./modules/bin_summary/main"

	// ANOTACION FUNCIONAL

include {METACERBERUS_CONTIGS} from "./modules/metacerberus/main"

	// PREDICCION TAXONOMICA EN CONTIGS

include {NT_BLASTN} from "./modules/nt_blastn/main"
include {BLOBTOOLS} from "./modules/blobtools/main"
include {SAMTOOLS_INDEX} from "./modules/samtools_index/main"
include {BLOBPLOT} from "./modules/blobplot/main"

	// REPORTES DE LECTURAS

include {READS_REPORT} from "./modules/reads_report/main"
include {TAX_REPORT_KRAKEN2} from "./modules/tax_report_kraken2/main"
include {TAX_REPORT_SOURMASH} from "./modules/tax_report_sourmash/main"

	// PREDICCION TAXONOMICA DE BINS

include {GTDB_TK} from "./modules/gtdb-tk/main"
include {GTDB_TK_COASSEMBLY} from "./modules/gtdb-tk/main"

	// CALCULO DE COVERTURA POR BIN

include {BEDTOOLS} from "./modules/bedtools/main"

	// PREDICCION DE ORFS EN CONTIGS Y BINS

include {PRODIGAL_BINS} from "./modules/prodigal/main"
include {PRODIGAL_CONTIGS} from "./modules/prodigal/main"

        // PREDICCION DE ARGS EN READS

include {KARGVA} from "./modules/kargva/main"
include {KARGA} from "./modules/karga/main"
include {ARGS_OAP} from "./modules/args_oap/main"
include {ARG_NORM_REPORT} from "./modules/arg_norm_report/main"

	// PREDICCION DE ARGS EN CONTIGS Y BINS

include {DEEPARG_BINS} from "./modules/deeparg/main"
include {DEEPARG_CONTIGS} from "./modules/deeparg/main"
include {ARG_CONTIG_LEVEL_REPORT} from "./modules/arg_contig_level_report/main"
include {ARG_FASTA_FORMATTER} from "./modules/arg_fasta_formatter/main"
include {CLUSTERING} from "./modules/clustering/main"
include {ARG_BLOBPLOT} from "./modules/arg_blobplot/main"

// Run the workflow

workflow{
    if ( params.help ) {

      help_info ='''
\u001B[33mUsage:\u001B[0m

\u001B[32mThe typical command for running the pipeline is as follows:\u001B[0m

nextflow run main.nf --input "path/to/samples_sheet" --output "path/to/output" -resume (recomended)

\u001B[33mMandatory arguments:\u001B[0m
   \u001B[32m--input\u001B[0m                        Input csv file with: samples names, path of all fastq files, and optionaly singletons.
                                  colnames required: "sample,r1,r2,s" if don't have singletons colname "s" can be empty

   \u001B[32m--output\u001B[0m                       Path to output dir 

\u001B[33mOptional arguments:\u001B[0m
   \u001B[32m--assembly_mode\u001B[0m                Mode of assembly, avaible options: "coassembly", "assembly", "none" 
                                  (default: assembly)
                                  coassembly: all samples are processing in one only data set
                                  assembly: all samples are processing individualy 

   \u001B[32m--taxonomic_profiler\u001B[0m           Software for taxonomic profiling, avaible options: "kraken2", "sourmash", "none"
                                  (default: sourmash)

   \u001B[32m--arg_clustering\u001B[0m               ARG gene prediction and clustering for horizontal gene transfer inference 
                                  (default: false)

   \u001B[32m--read_arg_prediction\u001B[0m          ARG and ARGV gene prediction at read level using KARGA and KARGVA 
                                  (default: false)

   \u001B[32m--contig_tax_and_arg\u001B[0m           Tax and ARG gene prediction at contig level using Blobtools and Deeparg
                                  (default: false)

   \u001B[32m--contig_level_metacerberus\u001B[0m    Contig level functional annotation with metacebeus 
                                  (default: false)

   \u001B[32m--help\u001B[0m                         Print this usage statement.

\u001B[33mAdditionally, all options can be modified in nextflow.config file\u001B[0m

   '''
    println(help_info)

    exit0
}

    switch ( params.taxonomic_profiler ) {

       case "kraken2":

       if ( params.custom_kraken_db?.trim() ) {
             ch_kraken_ref_db_formated = Channel.from(path(params.custom_kraken_db))
       } else {
             ch_kraken_ref_db = Channel.fromList(params.kraken_ref_db[params.kraken2_db]["file"])
             ch_kraken_ref_db_formated = FORMAT_KRAKEN_DB(ch_kraken_ref_db)
       }

       break;

       case "sourmash":

       if ( params.custom_sourmash_db?.trim() ) {
             ch_sourmash_db = Channel.fromList(params.custom_sourmash_db).map { file(it) }
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
       if ( params.custom_karga_db?.trim() ) {
             ch_karga_db = Channel.from(file(params.custom_karga_db))
       } else {
             ch_karga_db = Channel.fromList(params.karga_ref_db[params.karga_db]["file"]).map { file(it) }
       }
    

    if ( params.custom_kargva_db?.trim() ) {
             ch_kargva_db = Channel.from(file(params.custom_kargva_db))
       } else {
             ch_kargva_db = Channel.fromList(params.kargva_ref_db[params.kargva_db]["file"]).map { file(it) }
       }
    }

    if ( params.quality_control ) {

       if ( params.custom_phiX_index?.trim() ) {
            ch_bowtie_phiX_index_formated = Channel.fromPath(params.custom_phiX_index)
       } else {
            ch_bowtie_phiX_index = Channel.fromList(params.bowtie_ref_genomes_for_build[params.phiX_index]["file"]).map { file(it) }
            ch_bowtie_phiX_index_formated = BUILD_PHIX_BOWTIE2_INDEX(ch_bowtie_phiX_index)
       }

       if ( params.custom_bowtie_host_index?.trim() ) {
            ch_bowtie_host_index_formated = Channel.from(path(params.custom_bowtie_host_index))
       } else {
            ch_bowtie_host_index = Channel.fromList(params.bowtie_ref_host_index[params.host_db]["file"])
            ch_bowtie_host_index_formated = FORMAT_BOWTIE_HOST(ch_bowtie_host_index)
       }
    }

    if ( params.contig_tax_and_arg ) {

       if ( params.custom_deeparg_db?.trim() ) {
           ch_deeparg_db = Channel.fromPath(params.custom_deeparg_db)
       } else {
           ch_deeparg_db = DOWNLOAD_DEEPARG_DB()
       }

       if ( params.custom_blast_db?.trim() ) {
           ch_blast_db_formated = Channel.fromPath(params.custom_blast_db)
       } else {
           ch_blast_db =  Channel.fromList(params.blast_ref_db[params.blast_db]["url"])
           ch_blast_db_formated = FORMAT_NT_BLAST_DB(ch_blast_db)
       }

       if ( params.custom_taxdump_files?.trim() ) {
           ch_taxdump_files_formated = Channel.from(path(params.custom_taxdump_files))
       } else {
           ch_taxdump_files = Channel.fromList(params.taxonomy_files[params.taxdump_files]["url"])
           ch_taxdump_files_formated = FORMAT_TAXDUMP_FILES(ch_taxdump_files)
       }
    }
    
    if ( params.include_binning ) {

       if ( params.custom_gtdbtk_db?.trim() ) {
           ch_gtdbtk_db_formated = Channel.fromPath(params.custom_gtdbtk_db)
       } else {
           ch_gtdbtk_db = Channel.fromList(params.gtdbtk_ref_db[params.gtdbtk_db]["url"])
           ch_gtdbtk_db_formated = DOWNLOAD_GTDBTK_DB(ch_gtdbtk_db)
       }

       if ( params.custom_checkm2_db?.trim() ) {
           ch_checkm2_db_formated = Channel.fromPath(params.custom_checkm2_db)
       } else {
           ch_checkm2_db = Channel.fromList(params.checkm2_ref_db[params.checkm2_db]["url"])
           ch_checkm2_db_formated = FORMAT_CHECKM2_DB(ch_checkm2_db)
       }

    }

	// LEER EL ARCHIVO DE MUESTRAS Y GENERAR EL CANAL
    
    ch_reads = Channel
                   .from(file(params.input))
                   .splitCsv(header:true)
                   .map { row ->
                       def meta = [:]
                           meta.id = row.sample
                       def r1 = row.r1
                       def r2 = row.r2
                       if( row.s ) {
                            def s = row.s
                            return [meta, [r1, r2, s]]
                       }
                       else {
                            return [meta, [r1, r2]]
                       }
                   }

      if ( params.quality_control ) {

	        // LIMPIEZA DE LECTURAS POR CALIDAD

	    ch_fastp_reads = FASTP(ch_reads)

	        // EXTRAIGO Y FORMATEO REPORTES

	    ch_fastp_reads_report = QFILTER(ch_fastp_reads)

		// EXTRAIGO EL NUMERO DE LECTURAS

	    ch_fastp_reads_num = ch_fastp_reads_report.qfilter.map{ 
				       def after_reads = it[2].text
				       if( Integer.parseInt(after_reads) >= params.min_read_sample ) {
						def reads = it[1]
						def meta = it[0]
						return [meta, reads]
				       }
				       else {
						return 
				       }
				}

	        // MAPEO DE LECTURAS CONTRA BASES DE DATOS
		// SE INTRODUCE EL PATH DEL INDEX, EL INDEX BASENAME Y UN ALIAS DE LA BASE DE DATOS PARA LA SALIDA QUE NO CONTENGA PUNTOS

	    ch_phyx_clean_reads = BOWTIE2_PHYX(ch_fastp_reads_num.combine(ch_bowtie_phiX_index_formated), "phiX")
	    ch_host_clean_reads = BOWTIE2_HOST(ch_phyx_clean_reads.reads.combine(ch_bowtie_host_index_formated), "host")

	        // UNIFICACIÓN DE REPORTES

	    ch_reads_report = READS_REPORT(ch_phyx_clean_reads.report
	                                .concat(ch_host_clean_reads.report)
	                                .concat(ch_fastp_reads_report.reads_report)
	                                .collect(),"phiX host")
       } else {

            ch_host_clean_reads = COUNT_READS(ch_reads)

            ch_reads_report = READS_REPORT(ch_host_clean_reads.reads_report
                                        .collect(), "none")
       }

    switch ( params.taxonomic_profiler ) {
       case "kraken2":

            ch_k2_taxonomy = KRAKEN2(ch_host_clean_reads.reads.combine(ch_kraken_ref_db_formated), params.kraken_db_used)

            TAX_REPORT_KRAKEN2(ch_k2_taxonomy.report
                           .concat(ch_reads_report)
                           .collect())

        // ESTIMACION DE ABUNDANCIA, SE DEBE ESPECIFICAR EL PATH DE LA BASE DE DATOS,
        // UN ALIAS PARA LA SALIDA, EL TAMAÑO MINIMO DE READ Y EL NIVEL TAXONOMICO
        // SI FALLA ESTA PARTE, SE DEBE CAMBIAR EL NIVEL TAXONOMICO

            ch_taxonomy_estimation = BRACKEN(ch_k2_taxonomy.kraken.combine(ch_kraken_ref_db_formated), params.kraken_db_used)

        // GENERAR ARCHIVO .BIOM

            ch_kraken_biom = KRAKEN_BIOM(ch_taxonomy_estimation.collect(), params.kraken_db_used)

        // GENERAR OBJETO PHYLOSEQ Y GRAFICOS DE ABUNDANCIA

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

	// PREDICCION DE ARGS EN READS

    if ( params.read_arg_prediction) {

          ch_args_oap = ARGS_OAP(ch_host_clean_reads.reads)
          ch_argv_prediction = KARGVA(ch_host_clean_reads.reads.combine(ch_kargva_db))
          ch_arg_prediction = KARGA(ch_argv_prediction.kargva_reads.combine(ch_karga_db))
          ARG_NORM_REPORT(ch_arg_prediction.concat(ch_argv_prediction.kargva_reports).concat(ch_args_oap).collect())

    }

    switch ( params.assembly_mode ) {
       case "coassembly": 

        // ENSAMBLE DE LECTURAS, FILTRADO DE CONTIGS Y CALCULO DE PROFUNDIDAD

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

        // BINNING Y REFINAMIENTO DE BINS

              ch_metabat2 = METABAT2_COASSEMBLY(ch_filtered_contigs_co.contigs.combine(ch_depth))
              ch_semibin = SEMIBIN_COASSEMBLY(ch_bowtie2_samtools.collect().combine(ch_filtered_contigs_co.contigs))
              ch_comebin = COMEBIN_COASSEMBLY(ch_bowtie2_samtools.collect().combine(ch_filtered_contigs_co.contigs))

              ch_all_bins = ch_metabat2.combine(ch_semibin).combine(ch_comebin)
              ch_metawrap_co = METAWRAP_COASSEMBLY(ch_all_bins)

        // ESTIMACIÓN DE CALIDAD, ASIGNACIÓN TAXONOMICA Y GENERACIÓN DE REPORTES

              ch_checkm = CHECKM2_COASSEMBLY(ch_all_bins.combine(ch_metawrap_co).combine(ch_checkm2_db_formated))
              ch_gtdb_tk = GTDB_TK_COASSEMBLY(ch_metawrap_co.combine(ch_gtdbtk_db_formated))
              ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(ch_host_clean_reads.reads.combine(ch_metawrap_co))
              ch_bin_cov = BEDTOOLS(ch_bin_depth)

              BIN_SUMMARY(ch_bin_cov.collect().combine(ch_gtdb_tk.report.collect()).combine(ch_checkm.metawrap_report.collect()))
         }

         break;
         
         case "assembly":

        // ENSAMBLE DE LECTURAS, FILTRADO DE CONTIGS Y CALCULO DE PROFUNDIDAD

         ch_assembly = MEGAHIT(ch_host_clean_reads.reads)
         ch_filtered_contigs = BBMAP(ch_assembly)

        // BINNING Y REFINAMIENTO DE BINS

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

        // ESTIMACIÓN DE CALIDAD Y ASIGNACIÓN TAXONOMICA

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

	// PREDICCION ARG EN BINS Y CLUSTERING

    if ( params.arg_bin_clustering) {

        ch_raw_orfs = PRODIGAL_BINS(ch_metawrap)
        ch_deeparg = DEEPARG_BINS(ch_raw_orfs.combine(ch_deeparg_db))     
        ch_arg_fasta = ARG_FASTA_FORMATTER(ch_raw_orfs.join(ch_deeparg))        
        ch_clusters = CLUSTERING(ch_arg_fasta.collect())
    }

}

workflow.onComplete = {
    // any workflow property can be used here
}
