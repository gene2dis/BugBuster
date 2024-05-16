// DSL-2 enabled
nextflow.enable.dsl = 2

/*
IMPORT MODULES
*/

text='''
\u001B[36m    _________________   \u001B[0m    \u001B[33m__  __\u001B[0m \u001B[36m _                _     _       _\u001B[0m
\u001B[36m  / -.--.------.--.-- \\\u001B[0m    \u001B[33m|  \\/  |\u001B[0m\u001B[36m(_) ___ _ __ ___ | |__ (_) __ _| |\u001B[0m
\u001B[36m /-/  o--o\u001B[0m     \u001B[33mo--o\u001B[0m  \u001B[36m\\-\\\u001B[0m   \u001B[33m| |\\/| |\u001B[0m\u001B[36m| |/ __| °__/ _ \\|  _ \\| |/ _° | |\u001B[0m
\u001B[36m/-/  o--o\u001B[0m       \u001B[33mo--o\u001B[0m  \u001B[36m\\-\\\u001B[0m  \u001B[33m| |  | |\u001B[0m\u001B[36m| | (__| | | (_) | |_) | | (_| | |\u001B[0m
\u001B[36m|-|   o--o\u001B[0m     \u001B[33mo--o\u001B[0m   \u001B[36m|-|\u001B[0m  \u001B[33m|_|  |_|\u001B[0m\u001B[36m|_|\\___|_|  \\___/|_.__/|_|\\____|_|\u001B[0m
\u001B[36m|-|    o--o\u001B[0m   \u001B[33mo--o\u001B[0m    \u001B[36m|-|\u001B[0m   \u001B[33m____\u001B[0m\u001B[36m      °  _    °     \u001B[0m\u001B[33m____ \u001B[0m\u001B[36m    °  _         °\u001B[0m
\u001B[36m|-|\u001B[0m       \u001B[33mo---o\u001B[0m       \u001B[36m|-|\u001B[0m  \u001B[33m|  _ \\\u001B[0m\u001B[36m ° __ _| |_ __ _ °\u001B[0m\u001B[33m/ ___| \u001B[0m\u001B[36m  ___(_) ___°_ __  °___ ___\u001B[0m
\u001B[36m|-|\u001B[0m     \u001B[33mo--o\u001B[0m \u001B[36mo--o\u001B[0m     \u001B[36m|-|\u001B[0m  \u001B[33m| | | | \u001B[0m\u001B[36m/ _° | __/ _° | \u001B[0m\u001B[33m\\___ \\ \u001B[0m\u001B[36m / __| |/ _ \\ °_ \\ / __/ _ \\\u001B[0m
\u001B[36m|-|\u001B[0m    \u001B[33mo--o\u001B[0m   \u001B[36mo--o\u001B[0m    \u001B[36m|-|\u001B[0m  \u001B[33m| |_| |\u001B[0m\u001B[36m| (_| | || (_| |  \u001B[0m\u001B[33m___) |\u001B[0m\u001B[36m| (__| |  __/ | | | (_|  __/\u001B[0m
\u001B[36m|-|\u001B[0m   \u001B[33mo--o\u001B[0m     \u001B[36mo--o\u001B[0m   \u001B[36m|-|\u001B[0m  \u001B[33m|____/ \u001B[0m\u001B[36m \\____|\\__\\____| \u001B[0m\u001B[33m|____/ \u001B[0m\u001B[36m \\___|_|\\___|_| |_|\\___\\___|\u001B[0m
'''

println text 

include {FASTP} from "./modules/fastp/main"
include {QFILTER} from "./modules/qfilter/main"

include {BOWTIE2 as BOWTIE2_HUMAN} from "./modules/bowtie2/main"
include {BOWTIE2 as BOWTIE2_PHYX} from "./modules/bowtie2/main"

include {KRAKEN2 as KRAKEN2_GTDB} from "./modules/kraken2/main"
// include {KRAKEN2 as KRAKEN2_VIRAL} from "./modules/kraken2/main"
// include {KRAKEN2 as KRAKEN2_NT} from "./modules/kraken2/main"

include {BRACKEN as BRACKEN_GTDB} from "./modules/bracken/main"
// include {BRACKEN as BRACKEN_VIRAL} from "./modules/bracken/main"
// include {BRACKEN as BRACKEN_NT} from "./modules/bracken/main"

include {KRAKEN_BIOM as KRAKEN_BIOM_GTDB} from "./modules/kraken_biom/main"
// include {KRAKEN_BIOM as KRAKEN_BIOM_VIRAL} from "./modules/kraken_biom/main"
// include {KRAKEN_BIOM as KRAKEN_BIOM_NT} from "./modules/kraken_biom/main"

include {KRAKEN_TO_PHYLOSEQ} from "./modules/kraken_to_phyloseq/main"

include {TAX_REPORT} from "./modules/tax_report/main"

include {MEGAHIT} from "./modules/megahit/main_sin"
include {MEGAHIT as MEGAHIT_COASSEMBLY} from "./modules/megahit/main_co"

include {METACERBERUS_READS} from "./modules/metacerberus/main"
include {METACERBERUS_CONTIGS} from "./modules/metacerberus/main"

include {BBMAP} from "./modules/bbmap/main"
include {BBMAP as BBMAP_COASSEMBLY} from "./modules/bbmap/main_co"

include {BOWTIE2_SAMTOOLS} from "./modules/bowtie2_samtools/main"
include {BOWTIE2_SAMTOOLS as BOWTIE2_SAMTOOLS_COASSEMBLY} from "./modules/bowtie2_samtools/main_co"
include {BOWTIE2_SAMTOOLS as BOWTIE2_SAMTOOLS_DEPTH} from "./modules/bowtie2_samtools/main_depth"

include {CALCULATE_DEPTH} from "./modules/calculate_depth/main"
include {CALCULATE_DEPTH as CALCULATE_DEPTH_COASSEMBLY} from "./modules/calculate_depth/main_co"

include {SEMIBIN} from "./modules/semibin/main"
include {SEMIBIN as SEMIBIN_COASSEMBLY} from "./modules/semibin/main_co"

include {METABAT2} from "./modules/metabat2/main"
include {METABAT2 as METABAT2_COASSEMBLY} from "./modules/metabat2/main_co"

include {AUTOMETA} from "./modules/autometa/main"
include {AUTOMETA as AUTOMETA_COASSEMBLY} from "./modules/autometa/main_co"

include {METAWRAP} from "./modules/metawrap/main"
include {METAWRAP as METAWRAP_COASSEMBLY} from "./modules/metawrap/main_co"

include {CHECKM2} from "./modules/checkm2/main"
include {CHECKM2 as CHECKM2_COASSEMBLY} from "./modules/checkm2/main_co"

include {BIN_QUALITY_REPORT} from "./modules/bin_quality_report/main"
include {BIN_TAX_REPORT} from "./modules/bin_tax_report/main"
include {BIN_SUMMARY} from "./modules/bin_summary/main"

include {KARGVA} from "./modules/kargva/main"
include {KARGA} from "./modules/karga/main"
include {READS_REPORT} from "./modules/reads_report/main"

include {GTDB_TK} from "./modules/gtdb-tk/main"
include {GTDB_TK as GTDB_TK_COASSEMBLY} from "./modules/gtdb-tk/main_co"

include {BEDTOOLS} from "./modules/bedtools/main"

include {PRODIGAL} from "./modules/prodigal/main"
include {MMSEQS2} from "./modules/mmseqs2/main"
include {BEST_HITS_TO_FASTA} from "./modules/best_hits_to_fasta/main"
include {CLUSTERING} from "./modules/clustering/main"

// include {DARKHORSE} from "./modules/darkhorse/main"


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
   \u001B[32m--assembly_mode\u001B[0m                Mode of assembly, avaible options: "coassembly", "assembly" 
                                  (default: assembly)
                                  coassembly: all samples are processing in one only data set
                                  assembly: all samples are processing individualy 

   \u001B[32m--single_assembly\u001B[0m              If coassembly type is choosen, single assembly will aditionaly generate indivdual assembly for all samples 
                                  (default: false)

   \u001B[32m--arg_clustering\u001B[0m               ARG gene prediction and clustering for horizontal gene transfer inference 
                                  (default: false)

   \u001B[32m--read_arg_prediction\u001B[0m          ARG and ARGV gene prediction at read level using KARGA and KARGVA 
                                  (default: false)

   \u001B[32m--read_level_metacerberus\u001B[0m      Read level functional annotation with metacerberus 
                                  (default: false)

   \u001B[32m--contig_level_metacerberus\u001B[0m    Contig level functional annotation with metacebeus 
                                  (default: false)

   \u001B[32m--help\u001B[0m                         Print this usage statement.

\u001B[33mAdditionally, all options can be modified in nextflow.config file\u001B[0m

   '''
        println(help_info)
        exit 0
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
                       if( row.s != null ) {
                            def s = row.s
                            return [meta, [r1, r2, s]]
                       }
                       else {
                            return [meta, [r1, r2]]
                       }
                   }


        // LIMPIEZA DE LECTURAS POR CALIDAD

    ch_fastp_reads = FASTP(ch_reads)

        // EXTRAIGO Y FORMATEO REPORTES

    ch_fastp_reads_report = QFILTER(ch_fastp_reads)

	// EXTRAIGO EL NUMERO DE LECTURAS

    ch_fastp_reads_num = ch_fastp_reads_report.qfilter.map{ 
			       def after_reads = it[2].text
			       if( Integer.parseInt(after_reads[0]) >= 1000000 ) {
					def reads = it[1]
					def meta = it[0]
					return [meta, reads]
			       }
			       else {
					def reads = it[1]
					def meta = it[0]
					return [meta, reads]
			       }
			  }

        // MAPEO DE LECTURAS CONTRA BASES DE DATOS

	// SE INTRODUCE EL PATH DEL INDEX, EL INDEX BASENAME Y UN ALIAS DE LA BASE DE DATOS PARA LA SALIDA QUE NO CONTENGA PUNTOS

    ch_phyx_clean_reads = BOWTIE2_PHYX(ch_fastp_reads_num, params.phyX_db, "phiX_refseq", "phiX")
    ch_human_clean_reads = BOWTIE2_HUMAN(ch_phyx_clean_reads.reads, params.human_db, "chm13.draft_v1.0_plusY", "human")

        // UNIFICACIÓN DE REPORTES

    ch_reads_report = READS_REPORT(ch_phyx_clean_reads.report
				    .concat(ch_human_clean_reads.report)
				    .concat(ch_fastp_reads_report.reads_report)
				    .collect())

	// PREDICCIÓN TAXONOMICA

    ch_taxonomy_gtdb = KRAKEN2_GTDB(ch_human_clean_reads.reads, params.k2_gtdb_db, "gtdb")
//    ch_taxonomy_viral = KRAKEN2_VIRAL(ch_human_clean_reads.reads, params.k2_viral_db, "viral")
//    ch_taxonomy_nt = KRAKEN2_NT(ch_human_clean_reads.reads, params.k2_nt_db, "nt")


	// UNIFICACION DE REPORTES

    TAX_REPORT(ch_taxonomy_gtdb.report
			   .concat(ch_reads_report)
//			   .concat(ch_taxonomy_viral.report)
//			   .concat(ch_taxonomy_nt.report)
			   .collect())


	// ESTIMACION DE ABUNDANCIA, SE DEBE ESPECIFICAR EL PATH DE LA BASE DE DATOS,
	// UN ALIAS PARA LA SALIDA, EL TAMAÑO MINIMO DE READ Y EL NIVEL TAXONOMICO
	// SI FALLA ESTA PARTE, SE DEBE CAMBIAR EL NIVEL TAXONOMICO

    ch_taxonomy_estimation_gtdb = BRACKEN_GTDB(ch_taxonomy_gtdb.kraken, params.k2_gtdb_db, "gtdb")
//    ch_taxonomy_estimation_viral = BRACKEN_VIRAL(ch_taxonomy_viral.kraken, params.k2_viral_db, "viral")
//    ch_taxonomy_estimation_nt = BRACKEN_NT(ch_taxonomy_nt.kraken, params.k2_nt_db, "nt")

	// GENERAR ARCHIVO .BIOM
    ch_kraken_biom = KRAKEN_BIOM_GTDB(ch_taxonomy_estimation_gtdb.collect(), "gtdb")
//    KRAKEN_BIOM_VIRAL(ch_taxonomy_estimation_viral.collect(), "viral")
//    KRAKEN_BIOM_NT(ch_taxonomy_estimation_nt.collect(), "nt")

	// GENERAR OBJETO PHYLOSEQ Y GRAFICOS DE ABUNDANCIA

    KRAKEN_TO_PHYLOSEQ(ch_kraken_biom)

    switch ( params.assembly_mode ) {
       case "coassembly": 

        // ENSAMBLE DE LECTURAS, FILTRADO DE CONTIGS Y CALCULO DE PROFUNDIDAD

         ch_assembly_co = MEGAHIT_COASSEMBLY(ch_human_clean_reads.reads_coassembly.collect())
         ch_filtered_contigs_co = BBMAP_COASSEMBLY(ch_assembly_co)
         ch_bowtie2_samtools = BOWTIE2_SAMTOOLS_COASSEMBLY(ch_human_clean_reads.reads.combine(ch_filtered_contigs_co.contigs))
         ch_depth = CALCULATE_DEPTH_COASSEMBLY(ch_bowtie2_samtools.collect())

        // BINNING Y REFINAMIENTO DE BINS

         ch_metabat2 = METABAT2_COASSEMBLY(ch_filtered_contigs_co.contigs.combine(ch_depth))
         ch_semibin = SEMIBIN_COASSEMBLY(ch_bowtie2_samtools.collect().combine(ch_filtered_contigs_co.contigs))
         ch_autometa = AUTOMETA_COASSEMBLY(ch_bowtie2_samtools.collect().combine(ch_filtered_contigs_co.contigs), params.ncbi_db)

         ch_all_bins = ch_metabat2.combine(ch_semibin).combine(ch_autometa)
         ch_metawrap_co = METAWRAP_COASSEMBLY(ch_all_bins, params.metawrap_db)

        // ESTIMACIÓN DE CALIDAD, ASIGNACIÓN TAXONOMICA Y GENERACIÓN DE REPORTES

         ch_checkm = CHECKM2_COASSEMBLY(ch_all_bins.combine(ch_metawrap_co), params.checkm_db)
         ch_gtdb_tk = GTDB_TK_COASSEMBLY(ch_metawrap_co, params.gtdbtk_db)
         ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(ch_human_clean_reads.reads.combine(ch_metawrap_co))
         ch_bin_cov = BEDTOOLS(ch_bin_depth)

         BIN_SUMMARY(ch_bin_cov.collect().combine(ch_gtdb_tk.report.collect()).combine(ch_checkm.metawrap_report.collect()))

         if ( params.single_assembly ) {

              ch_assembly = MEGAHIT(ch_human_clean_reads.reads)
              ch_filtered_contigs = BBMAP(ch_assembly)

                if ( params.contig_level_metacerberus) {
                     METACERBERUS_CONTIGS(ch_filtered_contigs.bbmap_contigs)
                }
         }
         break;
         
         case "assembly":

        // ENSAMBLE DE LECTURAS, FILTRADO DE CONTIGS Y CALCULO DE PROFUNDIDAD

         ch_assembly = MEGAHIT(ch_human_clean_reads.reads)
         ch_filtered_contigs = BBMAP(ch_assembly)
         ch_bowtie2_samtools = BOWTIE2_SAMTOOLS(ch_filtered_contigs.bbmap)
         ch_depth = CALCULATE_DEPTH(ch_bowtie2_samtools)

        // BINNING Y REFINAMIENTO DE BINS

         ch_metabat2 = METABAT2(ch_depth)
         ch_semibin = SEMIBIN(ch_bowtie2_samtools)
         ch_autometa = AUTOMETA(ch_bowtie2_samtools, params.ncbi_db)

         ch_all_bins = ch_metabat2.join(ch_semibin).join(ch_autometa)
         ch_metawrap = METAWRAP(ch_all_bins, params.metawrap_db)

         ch_filtered_contigs.bbmap.join(ch_metawrap).join(ch_human_clean_reads.reads).view()

        // ESTIMACIÓN DE CALIDAD Y ASIGNACIÓN TAXONOMICA

         ch_checkm = CHECKM2(ch_all_bins.join(ch_metawrap), params.checkm_db)
         ch_gtdb_tk = GTDB_TK(ch_metawrap, params.gtdbtk_db)   


         if ( params.contig_level_metacerberus) {
               METACERBERUS_CONTIGS(ch_filtered_contigs.bbmap_contigs)
         }  
         break;

         default:
         
         exit 1, "ERROR: Not valid assembly type, check nextflow.config file "
    }

	// UNIFICACIÓN DE REPORTES

    BIN_QUALITY_REPORT(ch_checkm.all_reports.collect())
    BIN_TAX_REPORT(ch_gtdb_tk.report.collect())

	// PREDICCION DE ORF

    if ( params.arg_clustering) {

        ch_raw_orfs = PRODIGAL(ch_metawrap)
        ch_mmseqs = MMSEQS2(ch_raw_orfs, params.arg_db)     
    
        ch_best_hits_faa = BEST_HITS_TO_FASTA(ch_mmseqs.collect())
        ch_clusters = CLUSTERING(ch_best_hits_faa)
    }

	// PREDICCIÓN DE ARGV y ARG A NIVEL DE READS

    if ( params.read_arg_prediction) {

        ch_argv_prediction = KARGVA(ch_human_clean_reads.reads) 
        ch_arg_prediction = KARGA(ch_argv_prediction, params.karga_db)
    }
}

workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Vamos bien mi king"
}
