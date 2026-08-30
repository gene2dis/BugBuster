/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AGGREGATE_FUNCTIONS Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Study-level aggregation of the contig functional branch (design doc
    Section 4.8, tasks T4+T5): joins featureCounts gene counts, eggNOG-mapper
    v3 annotations and Pyrodigal gene coordinates into the canonical tables
    (Sections 5.1-5.3), computes TPM (Section 6.1) and copies per genome
    equivalent (Section 6.2, from the optional MicrobeCensus AGS tables -
    samples without one fall back to TPM-only with an 'unavailable' row in
    ags_and_ge.tsv), and reports the annotated fraction per sample. All the
    logic lives in
    bin/aggregate_functions.py (independently tested by
    tests/bin/test_aggregate_functions.sh), including the version-aware
    eggNOG layout guard that fails loudly on schema drift (Section 4.2).

    Sample identity is derived from input filenames (<id>.featureCounts.txt
    etc.), which equal meta.id because every upstream module names outputs
    with task.ext.prefix ?: meta.id - do not set ext.prefix on the upstream
    functional modules without revisiting this.

    Input:
        counts: all samples' featureCounts tables
        annotations: eggNOG annotations (per sample, or one 'coassembly' file)
        gffs: Pyrodigal gene GFFs (per sample, or one 'coassembly' file)
        ags: MicrobeCensus <id>.ags.tsv tables; possibly a subset of samples
            or none at all (non-fatal failures, or --microbecensus false)
        eggnog_versions: versions.yml from EGGNOG_MAPPER_ANNOTATE (staged
            under a distinct name: this module writes its own versions.yml,
            and writing through a same-named input symlink would corrupt the
            upstream work directory - same hazard as TAXONOMY_REPORT's
            input_reads_report.csv)
        assembly_mode: 'assembly' or 'coassembly' (join topology + 5.2 column)

    Output: gene_annotations.tsv, gene_abundance.tsv, function_abundance.tsv,
        function_wide_<ontology>_{tpm,cpge}.tsv (ko, cog, ec, pfam, cazy),
        annotated_fraction.tsv, ags_and_ge.tsv, versions.yml
----------------------------------------------------------------------------------------
*/

process AGGREGATE_FUNCTIONS {
    tag "functions-${assembly_mode}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pruned:3b27e3935fded4e0' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib_pruned:195b3e3e5f741210' }"

    input:
    path(counts, stageAs: 'counts/*')
    path(annotations, stageAs: 'annotations/*')
    path(gffs, stageAs: 'gffs/*')
    path(ags, stageAs: 'ags/*')
    path(eggnog_versions, stageAs: 'eggnog_versions.yml')
    val(assembly_mode)

    output:
    path("gene_annotations.tsv")   , emit: gene_annotations
    path("gene_abundance.tsv")     , emit: gene_abundance
    path("function_abundance.tsv") , emit: function_abundance
    path("function_wide_*.tsv")    , emit: function_wide
    path("annotated_fraction.tsv") , emit: annotated_fraction
    path("ags_and_ge.tsv")         , emit: ags_summary
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ags_arg = ags ? "--ags ${ags}" : ''
    """
    aggregate_functions.py \\
        --assembly-mode ${assembly_mode} \\
        --counts ${counts} \\
        --annotations ${annotations} \\
        --gffs ${gffs} \\
        ${ags_arg} \\
        --eggnog-versions-yml ${eggnog_versions} \\
        --output-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch gene_annotations.tsv
    touch gene_abundance.tsv
    touch function_abundance.tsv
    touch function_wide_ko_tpm.tsv
    touch function_wide_cog_tpm.tsv
    touch function_wide_ec_tpm.tsv
    touch function_wide_pfam_tpm.tsv
    touch function_wide_cazy_tpm.tsv
    touch function_wide_ko_cpge.tsv
    touch function_wide_cog_cpge.tsv
    touch function_wide_ec_cpge.tsv
    touch function_wide_pfam_cpge.tsv
    touch function_wide_cazy_cpge.tsv
    touch annotated_fraction.tsv
    touch ags_and_ge.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
        pandas: 2.0
    END_VERSIONS
    """
}
