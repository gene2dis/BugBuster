/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AGGREGATE_READ_FUNCTIONS Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Study-level aggregation of the read-level functional branch (design doc
    Sections 4.6 and 5.3, tasks T8a/T8b): converts every sample's per-function
    read abundances into the canonical function abundance schema with
    source=reads, plus wide matrices, annotated fractions and a per-sample
    summary. Thin wrapper around bin/aggregate_read_functions.py (version-
    aware: fails on an unverified backend version or input layout).

    Reported SEPARATELY from the contig branch (Section 2): all outputs are
    read_*-prefixed and never merged into function_abundance.tsv.

    CPGE (RPK / genome equivalents, Section 6.2) needs the MicrobeCensus AGS
    tables, which are optional: samples whose MicrobeCensus failed (non-fatal)
    or runs with --microbecensus false get blank CPGE and an 'unavailable'
    cpge_status in read_sample_summary.tsv. The superfocus backend has no
    gene length: CPGE is 'not_applicable' and no _cpge matrices are written
    (design doc Q17).

    Input:  path(functions)   all samples' <id>.<backend>_functions.tsv
            path(summaries)   all samples' <id>.<backend>_summary.tsv
            path(unassigned)  all samples' <id>.woltka_unassigned.tsv
                              (woltka only; [] for superfocus)
            path(ags)         MicrobeCensus <id>.ags.tsv files (may be empty)
            path(backend_versions)  WOLTKA_CLASSIFY / SUPERFOCUS versions.yml
                              (staged under a distinct name: this process
                              emits its own)
            val(backend)      'woltka' | 'superfocus'
    Output: read_function_abundance.tsv, read_function_wide_*.tsv,
            read_annotated_fraction.tsv, read_sample_summary.tsv, versions.yml
----------------------------------------------------------------------------------------
*/

process AGGREGATE_READ_FUNCTIONS {
    tag "read-functions-${backend}"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pruned:3b27e3935fded4e0' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pruned:195b3e3e5f741210' }"

    input:
    path(functions, stageAs: 'functions/*')
    path(summaries, stageAs: 'summaries/*')
    path(unassigned, stageAs: 'unassigned/*')
    path(ags, stageAs: 'ags/*')
    path(backend_versions, stageAs: 'backend_versions.yml')
    val(backend)

    output:
    path("read_function_abundance.tsv"), emit: function_abundance
    path("read_function_wide_*.tsv")   , emit: function_wide
    path("read_annotated_fraction.tsv"), emit: annotated_fraction
    path("read_sample_summary.tsv")    , emit: sample_summary
    path("versions.yml")               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ags_arg = ags ? "--ags ${ags}" : ''
    def unassigned_arg = unassigned ? "--unassigned ${unassigned}" : ''
    """
    aggregate_read_functions.py \\
        --backend ${backend} \\
        --functions ${functions} \\
        --summaries ${summaries} \\
        ${unassigned_arg} \\
        ${ags_arg} \\
        --versions-yml ${backend_versions} \\
        --output-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def wide_files = backend == 'superfocus' ?
        ['seed_level1', 'seed_level2', 'seed_level3'].collect { o -> "read_function_wide_${o}_native.tsv" } :
        ['ko', 'ec', 'cog', 'pfam', 'metacyc'].collect { o -> "read_function_wide_${o}_native.tsv read_function_wide_${o}_cpge.tsv" }
    """
    touch read_function_abundance.tsv
    touch ${wide_files.join(' ')}
    touch read_annotated_fraction.tsv
    touch read_sample_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
        pandas: 2.0
    END_VERSIONS
    """
}
