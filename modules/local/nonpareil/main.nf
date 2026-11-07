/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NONPAREIL Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sequencing coverage and diversity of one sample (roadmap Phase 17, N1): Nonpareil
    k-mer mode on R1 of the clean reads, then the Nonpareil R package's curve fit
    (bin/nonpareil_fit.R) in the same image.

    R1 only: Nonpareil estimates coverage from the redundancy of comparable reads, and
    R2 adds no fragments of its own; the singleton file is not used either. The
    estimates are per sample, so they hold in any batch.

    k-mer mode has no seed (nonpareil's -r applies to alignment mode only), so the curve
    and its values change slightly from run to run.

    Failure is non-fatal: the withName block in config/modules.config retries OOM-class
    exits, then ignores the task, and the sample has no Nonpareil values. A curve the
    package rejects (too few reads) or a model that does not converge is not a failure:
    bin/nonpareil_fit.R reports the coverage with blank Nd and LR*.

    Input:  tuple val(meta), path(reads)  clean reads (R1/R2[/S]); R1 is reads[0]
    Output: tuple val(meta), path("*_nonpareil.tsv")      , emit: values
            tuple val(meta), path("*_nonpareil_curve.tsv"), emit: curve
            tuple val(meta), path("*.{npo,npl}")          , emit: raw
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process NONPAREIL {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/nonpareil:3.5.5--r44h077b44d_2' :
        'quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_nonpareil.tsv")      , emit: values
    tuple val(meta), path("${meta.id}_nonpareil_curve.tsv"), emit: curve
    tuple val(meta), path("${meta.id}.{npo,npl}")          , emit: raw
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    def r1 = reads instanceof List ? reads[0] : reads
    // -R: memory for the reads, in MiB (Nonpareil's default 1024 caps it otherwise); 90 % of the task's
    def ram_mib = task.memory ? Math.min((task.memory.toMega() * 0.9) as long, 4194303L) : 1024
    """
    nonpareil \\
        -s ${r1} \\
        -T kmer \\
        -f fastq \\
        -b ${prefix} \\
        -t ${task.cpus} \\
        -R ${ram_mib} \\
        ${args}

    nonpareil_fit.R ${prefix}.npo ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(nonpareil -V 2>&1 | sed -n 's/^Nonpareil v//p')
        r-nonpareil: \$(Rscript -e 'cat(as.character(packageVersion("Nonpareil")))')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    printf 'sample\\treads\\tread_length\\tal\\tkappa\\tcoverage\\tdiversity\\tlr_star\\tdepth_95pct_pairs\\tmodel_r\\thas_model\\tconsistent\\twarning\\n' > ${prefix}_nonpareil.tsv
    printf 'sample\\tkind\\tx_adj\\tx_pairs\\tcoverage\\tcoverage_p25\\tcoverage_p75\\n' > ${prefix}_nonpareil_curve.tsv
    touch ${prefix}.npo ${prefix}.npl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: 3.5.5
        r-nonpareil: 3.5.3
    END_VERSIONS
    """
}
