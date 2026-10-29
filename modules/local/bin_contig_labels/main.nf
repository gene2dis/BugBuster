process BIN_CONTIG_LABELS {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*_contig_bins.tsv"), emit: labels
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bin_contig_labels.py ${meta.id} ${bins}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'contig\\tbin\\n' > ${meta.id}_contig_bins.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.12
    END_VERSIONS
    """
}
