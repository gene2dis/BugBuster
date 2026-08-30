/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONAL ANNOTATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Contig-level functional annotation (design doc Section 3.4): the
    eggNOG-mapper branch (T2, two-stage DIAMOND search + orthology-transfer
    annotation of the shared Pyrodigal proteins), gene quantification (T3,
    featureCounts over the Pyrodigal gene coordinates), and study-level
    aggregation into the canonical abundance tables (T4 TPM, T5 CPGE via the
    MicrobeCensus AGS channel — optional: samples whose MicrobeCensus failed
    or with --microbecensus false contribute no ags.tsv and fall back to
    TPM-only). dbCAN and the MAG branch are added by later tasks (T6, T7).

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
include { AGGREGATE_FUNCTIONS    } from '../../modules/local/aggregate_functions/main'

workflow FUNCTIONAL_ANNOTATION {

    take:
    ch_proteins     // channel: [ meta, proteins.faa(.gz) ] from the shared gene-calling step
    ch_eggnog_db    // channel: eggNOG 7 data dir (emapper-3.0 layout, see config/databases.config)
    ch_gene_gff     // channel: [ meta, genes.gff(.gz) ] from the shared gene-calling step
    ch_counting_bam // channel: [ meta, bam ] per-sample reads-vs-contigs BAMs (ASSEMBLY.out.counting_bam)
    ch_ags          // channel: [ meta, ags.tsv ] MicrobeCensus AGS tables; may be empty or miss samples (non-fatal failures)

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

    // Study-level aggregation (design doc Section 4.8): all samples' counts,
    // the annotations and GFFs (per-sample or single coassembly elements),
    // the eggNOG versions.yml for the version-aware layout guard, and the
    // MicrobeCensus AGS tables (ifEmpty([]) because the channel legitimately
    // emits nothing when MicrobeCensus is off or every sample's run failed)
    AGGREGATE_FUNCTIONS(
        FEATURECOUNTS_GENES.out.counts.map { _meta, counts -> counts }.collect(),
        EGGNOG_MAPPER_ANNOTATE.out.annotations.map { _meta, annotations -> annotations }.collect(),
        ch_gene_gff.map { _meta, gff -> gff }.collect(),
        ch_ags.map { _meta, ags -> ags }.collect().ifEmpty([]),
        EGGNOG_MAPPER_ANNOTATE.out.versions.first(),
        params.assembly_mode
    )

    ch_versions = ch_versions.mix(
        EGGNOG_MAPPER_SEARCH.out.versions.first(),
        EGGNOG_MAPPER_ANNOTATE.out.versions.first(),
        FEATURECOUNTS_GENES.out.versions.first(),
        AGGREGATE_FUNCTIONS.out.versions
    )

    emit:
    annotations         = EGGNOG_MAPPER_ANNOTATE.out.annotations       // [ meta, *.emapper.annotations ]
    seed_orthologs      = EGGNOG_MAPPER_SEARCH.out.seed_orthologs      // [ meta, *.emapper.seed_orthologs ]
    gene_counts         = FEATURECOUNTS_GENES.out.counts               // [ meta, *.featureCounts.txt ] - per sample in both modes
    gene_counts_summary = FEATURECOUNTS_GENES.out.summary              // [ meta, *.featureCounts.txt.summary ]
    gene_annotations    = AGGREGATE_FUNCTIONS.out.gene_annotations     // gene_annotations.tsv (Section 5.1)
    gene_abundance      = AGGREGATE_FUNCTIONS.out.gene_abundance       // gene_abundance.tsv (Section 5.2)
    function_abundance  = AGGREGATE_FUNCTIONS.out.function_abundance   // function_abundance.tsv (Section 5.3)
    function_wide       = AGGREGATE_FUNCTIONS.out.function_wide        // function_wide_<ontology>_{tpm,cpge}.tsv
    annotated_fraction  = AGGREGATE_FUNCTIONS.out.annotated_fraction   // annotated_fraction.tsv
    ags_summary         = AGGREGATE_FUNCTIONS.out.ags_summary          // ags_and_ge.tsv
    versions            = ch_versions
}
