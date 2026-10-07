/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    READ FUNCTIONAL SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Optional read-level functional profiling of host-removed reads (design
    doc Sections 3.4 and 4.6). One backend per run, selected by
    --read_level_functional; independent of assembly (runs with
    --assembly_mode none) and of the contig/MAG functional branches.
    Read-level profiling recovers the unassembled fraction but over-predicts,
    so its tables are reported SEPARATELY (read_* files), never merged into
    the contig branch's function_abundance.tsv (Section 2).

    Backends:
      woltka (T8a, Section 4.6.2): WOLTKA_ALIGN (Bowtie2 vs WoLr2, SHOGUN
        multi-hit) -> WOLTKA_CLASSIFY (woltka classify + per-ORF term-set
        composition to ko/ec/cog/pfam/metacyc) -> AGGREGATE_READ_FUNCTIONS
        (canonical Section 5.3 schema, source=reads)

    CPGE uses the MicrobeCensus AGS channel (optional: samples whose
    MicrobeCensus failed, or --microbecensus false, keep blank CPGE).
----------------------------------------------------------------------------------------
*/

include { WOLTKA_ALIGN             } from '../../modules/local/woltka_align/main'
include { WOLTKA_CLASSIFY          } from '../../modules/local/woltka_classify/main'
include { AGGREGATE_READ_FUNCTIONS } from '../../modules/local/aggregate_read_functions/main'

workflow READ_FUNCTIONAL {

    take:
    ch_reads     // channel: [ meta, [ R1, R2(, Singleton) ] ] host-removed reads
    ch_woltka_db // channel: WoLr2-layout db dir; empty unless --read_level_functional woltka
    ch_ags       // channel: [ meta, ags.tsv ] MicrobeCensus AGS tables; may be empty or miss samples

    main:
    ch_versions = channel.empty()

    if ( params.read_level_functional == 'woltka' ) {
        WOLTKA_ALIGN(ch_reads.combine(ch_woltka_db))
        WOLTKA_CLASSIFY(WOLTKA_ALIGN.out.sam.combine(ch_woltka_db))

        AGGREGATE_READ_FUNCTIONS(
            WOLTKA_CLASSIFY.out.functions.map { _meta, f -> f }.collect(),
            WOLTKA_CLASSIFY.out.summary.map { _meta, f -> f }.collect(),
            WOLTKA_CLASSIFY.out.unassigned.map { _meta, f -> f }.collect(),
            ch_ags.map { _meta, ags -> ags }.collect().ifEmpty([]),
            WOLTKA_CLASSIFY.out.versions.first(),
            'woltka'
        )

        ch_sample_functions    = WOLTKA_CLASSIFY.out.functions
        ch_function_abundance  = AGGREGATE_READ_FUNCTIONS.out.function_abundance
        ch_function_wide       = AGGREGATE_READ_FUNCTIONS.out.function_wide
        ch_annotated_fraction  = AGGREGATE_READ_FUNCTIONS.out.annotated_fraction
        ch_sample_summary      = AGGREGATE_READ_FUNCTIONS.out.sample_summary

        ch_versions = ch_versions.mix(
            WOLTKA_ALIGN.out.versions.first(),
            WOLTKA_CLASSIFY.out.versions.first(),
            AGGREGATE_READ_FUNCTIONS.out.versions
        )
    } else {
        ch_sample_functions    = channel.empty()
        ch_function_abundance  = channel.empty()
        ch_function_wide       = channel.empty()
        ch_annotated_fraction  = channel.empty()
        ch_sample_summary      = channel.empty()
    }

    emit:
    sample_functions   = ch_sample_functions    // [ meta, <id>.<backend>_functions.tsv ] per sample
    function_abundance = ch_function_abundance  // read_function_abundance.tsv (Section 5.3, source=reads)
    function_wide      = ch_function_wide       // read_function_wide_<ontology>_{native,cpge}.tsv
    annotated_fraction = ch_annotated_fraction  // read_annotated_fraction.tsv
    sample_summary     = ch_sample_summary      // read_sample_summary.tsv
    versions           = ch_versions
}
