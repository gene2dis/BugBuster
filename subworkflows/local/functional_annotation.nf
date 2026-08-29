/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONAL ANNOTATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Contig-level functional annotation (design doc Section 3.4): the
    eggNOG-mapper branch (T2, two-stage DIAMOND search + orthology-transfer
    annotation of the shared Pyrodigal proteins) and gene quantification (T3,
    featureCounts over the Pyrodigal gene coordinates). Aggregation, dbCAN and
    the MAG branch are added by later tasks (T4-T7).

    Assembly-mode handling (design doc Section 3.5): annotation works
    identically in both modes (per-sample elements under 'assembly', a single
    'coassembly' element under 'coassembly'). Quantification is mode-aware:
    under 'assembly' each sample's BAM joins its own gene GFF; under
    'coassembly' every per-sample BAM is counted against the single shared
    co-assembly GFF, giving per-sample counts over a shared gene set.
----------------------------------------------------------------------------------------
*/

include { EGGNOG_MAPPER_SEARCH   } from '../../modules/local/eggnog_mapper_search/main'
include { EGGNOG_MAPPER_ANNOTATE } from '../../modules/local/eggnog_mapper_annotate/main'
include { FEATURECOUNTS_GENES    } from '../../modules/local/featurecounts_genes/main'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteins     // channel: [ meta, proteins.faa(.gz) ] from the shared gene-calling step
    ch_eggnog_db    // channel: eggNOG 7 data dir (emapper-3.0 layout, see config/databases.config)
    ch_gene_gff     // channel: [ meta, genes.gff(.gz) ] from the shared gene-calling step
    ch_counting_bam // channel: [ meta, bam ] per-sample reads-vs-contigs BAMs (ASSEMBLY.out.counting_bam)

    main:
    ch_versions = channel.empty()

    EGGNOG_MAPPER_SEARCH(ch_proteins.combine(ch_eggnog_db))
    EGGNOG_MAPPER_ANNOTATE(EGGNOG_MAPPER_SEARCH.out.seed_orthologs.combine(ch_eggnog_db))

    // Pair each counting BAM with its gene GFF (design doc Section 4.4):
    // under coassembly the GFF channel is a single 'coassembly' element shared
    // by every per-sample BAM, so a meta join cannot apply
    if ( params.assembly_mode == 'coassembly' ) {
        ch_fc_input = ch_counting_bam.combine(ch_gene_gff.map { _meta, gff -> gff })
    } else {
        ch_fc_input = ch_counting_bam.join(ch_gene_gff)
    }

    FEATURECOUNTS_GENES(ch_fc_input)

    ch_versions = ch_versions.mix(
        EGGNOG_MAPPER_SEARCH.out.versions.first(),
        EGGNOG_MAPPER_ANNOTATE.out.versions.first(),
        FEATURECOUNTS_GENES.out.versions.first()
    )

    emit:
    annotations         = EGGNOG_MAPPER_ANNOTATE.out.annotations   // [ meta, *.emapper.annotations ]
    seed_orthologs      = EGGNOG_MAPPER_SEARCH.out.seed_orthologs  // [ meta, *.emapper.seed_orthologs ]
    gene_counts         = FEATURECOUNTS_GENES.out.counts           // [ meta, *.featureCounts.txt ] - per sample in both modes
    gene_counts_summary = FEATURECOUNTS_GENES.out.summary          // [ meta, *.featureCounts.txt.summary ]
    versions            = ch_versions
}
