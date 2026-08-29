/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONAL ANNOTATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Contig-level functional annotation (design doc Section 3.4). Currently the
    eggNOG-mapper branch (T2): two-stage DIAMOND search + orthology-transfer
    annotation of the shared Pyrodigal proteins. Gene quantification,
    aggregation, dbCAN and the MAG branch are added by later tasks (T3-T7).

    Works identically under both assembly modes: per-sample proteins under
    'assembly', a single 'coassembly' element under 'coassembly'.
----------------------------------------------------------------------------------------
*/

include { EGGNOG_MAPPER_SEARCH   } from '../../modules/local/eggnog_mapper_search/main'
include { EGGNOG_MAPPER_ANNOTATE } from '../../modules/local/eggnog_mapper_annotate/main'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteins   // channel: [ meta, proteins.faa(.gz) ] from the shared gene-calling step
    ch_eggnog_db  // channel: eggNOG 7 data dir (emapper-3.0 layout, see config/databases.config)

    main:
    ch_versions = channel.empty()

    EGGNOG_MAPPER_SEARCH(ch_proteins.combine(ch_eggnog_db))
    EGGNOG_MAPPER_ANNOTATE(EGGNOG_MAPPER_SEARCH.out.seed_orthologs.combine(ch_eggnog_db))

    ch_versions = ch_versions.mix(
        EGGNOG_MAPPER_SEARCH.out.versions.first(),
        EGGNOG_MAPPER_ANNOTATE.out.versions.first()
    )

    emit:
    annotations    = EGGNOG_MAPPER_ANNOTATE.out.annotations   // [ meta, *.emapper.annotations ]
    seed_orthologs = EGGNOG_MAPPER_SEARCH.out.seed_orthologs  // [ meta, *.emapper.seed_orthologs ]
    versions       = ch_versions
}
