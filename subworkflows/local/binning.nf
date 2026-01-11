/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BINNING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic binning using multiple tools and refinement with MetaWRAP
----------------------------------------------------------------------------------------
*/

include { CALCULATE_DEPTH            } from '../../modules/local/calculate_depth/main'
include { CALCULATE_DEPTH_COASSEMBLY } from '../../modules/local/calculate_depth/main'
include { METABAT2_METABAT2 as METABAT2 } from '../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_COASSEMBLY        } from '../../modules/local/metabat2/main'
include { SEMIBIN                    } from '../../modules/local/semibin/main'
include { SEMIBIN_COASSEMBLY         } from '../../modules/local/semibin/main'
include { COMEBIN                    } from '../../modules/local/comebin/main'
include { COMEBIN_COASSEMBLY         } from '../../modules/local/comebin/main'
include { METAWRAP                   } from '../../modules/local/metawrap/main'
include { METAWRAP_COASSEMBLY        } from '../../modules/local/metawrap/main'
include { CHECKM2                    } from '../../modules/local/checkm2/main'
include { CHECKM2_COASSEMBLY         } from '../../modules/local/checkm2/main'
include { GTDB_TK                    } from '../../modules/local/gtdb-tk/main'
include { GTDB_TK_COASSEMBLY         } from '../../modules/local/gtdb-tk/main'
include { BOWTIE2_SAMTOOLS_DEPTH     } from '../../modules/local/bowtie2_samtools/main'
include { BEDTOOLS                   } from '../../modules/local/bedtools/main'
include { BIN_QUALITY_REPORT         } from '../../modules/local/bin_quality_report/main'
include { BIN_TAX_REPORT             } from '../../modules/local/bin_tax_report/main'
include { BIN_SUMMARY                } from '../../modules/local/bin_summary/main'

workflow BINNING {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    contigs         // channel: [ val(meta), path(reads), path(contigs) ] or path(contigs)
    bam             // channel: path(bam) or [ val(meta), path(bam) ]
    checkm2_db      // channel: path(checkm2_db)
    gtdbtk_db       // channel: path(gtdbtk_db)

    main:
    ch_versions = Channel.empty()
    ch_refined_bins = Channel.empty()

    //
    // Per-sample binning mode
    //
    if ( params.assembly_mode == "assembly" ) {
        // Calculate depth
        ch_depth = CALCULATE_DEPTH(bam)
        
        //
        // Run nf-core METABAT2
        // nf-core METABAT2_METABAT2 signature:
        //   input:  tuple val(meta), path(fasta), path(depth)
        //   output: tuple val(meta), path("*.fa.gz"), emit: fasta
        //
        ch_metabat2_input = ch_depth.map { meta, contigs_file, depth_file ->
            [ meta, contigs_file, depth_file ]
        }
        
        METABAT2(ch_metabat2_input)
        
        // Collect versions
        ch_versions = ch_versions.mix(METABAT2.out.versions.first())
        
        // nf-core outputs individual gzipped bins; collect them for downstream
        // Group bins by sample for compatibility with local binners
        ch_metabat2 = METABAT2.out.fasta
        
        // bam channel already has the correct structure: [ meta, contigs, bam ]
        ch_semibin = SEMIBIN(bam)
        ch_comebin = COMEBIN(bam)

        // Combine bins for refinement
        ch_all_bins = ch_metabat2.join(ch_semibin).join(ch_comebin)
        
        // Refine bins with MetaWRAP
        ch_metawrap = METAWRAP(ch_all_bins)
        ch_refined_bins = ch_metawrap

        // Quality assessment with CheckM2
        ch_checkm = CHECKM2(
            ch_all_bins.join(ch_metawrap).combine(checkm2_db)
        )

        // Taxonomic classification with GTDB-TK
        ch_gtdb_tk = GTDB_TK(ch_metawrap.combine(gtdbtk_db))

        // Generate reports
        BIN_QUALITY_REPORT(ch_checkm.all_reports.collect())
        BIN_TAX_REPORT(ch_gtdb_tk.report.collect())
    }

    //
    // Co-assembly binning mode
    //
    if ( params.assembly_mode == "coassembly" ) {
        // Calculate depth from all BAMs
        ch_depth_co = CALCULATE_DEPTH_COASSEMBLY(bam.collect())
        
        // Run multiple binners
        ch_metabat2_co = METABAT2_COASSEMBLY(
            contigs.combine(ch_depth_co)
        )
        ch_semibin_co = SEMIBIN_COASSEMBLY(
            bam.collect().combine(contigs)
        )
        ch_comebin_co = COMEBIN_COASSEMBLY(
            bam.collect().combine(contigs)
        )

        // Combine and refine bins
        ch_all_bins_co = ch_metabat2_co.combine(ch_semibin_co).combine(ch_comebin_co)
        ch_metawrap_co = METAWRAP_COASSEMBLY(ch_all_bins_co)
        ch_refined_bins = ch_metawrap_co

        // Quality assessment
        ch_checkm_co = CHECKM2_COASSEMBLY(
            ch_all_bins_co.combine(ch_metawrap_co).combine(checkm2_db)
        )

        // Taxonomic classification
        ch_gtdb_tk_co = GTDB_TK_COASSEMBLY(ch_metawrap_co.combine(gtdbtk_db))

        // Calculate bin coverage
        ch_bin_depth = BOWTIE2_SAMTOOLS_DEPTH(reads.combine(ch_metawrap_co))
        ch_bin_cov = BEDTOOLS(ch_bin_depth)

        // Generate summary report
        BIN_SUMMARY(
            ch_bin_cov.collect()
                .combine(ch_gtdb_tk_co.report.collect())
                .combine(ch_checkm_co.metawrap_report.collect())
        )
    }

    emit:
    refined_bins = ch_refined_bins    // channel: [ val(meta), path(bins) ] or path(bins)
    versions     = ch_versions        // channel: path(versions.yml)
}
