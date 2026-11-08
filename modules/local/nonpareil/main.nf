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

    k-mer mode stops ("Reads are required to have a minimum length of kmer size") on any
    read shorter than k, and fastp keeps reads down to 15 bp: R1 reads shorter than k
    (-k in ext.args, 24 by default) are dropped first, into an uncompressed copy in the
    task directory (the copy nonpareil makes of a .gz input anyway), removed at the end.
    Their number is the reads_shorter_than_k column (0.09 % of the R1 reads of a real
    10 M-pair sample).

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
    def kmer = (args =~ /-k\s+(\d+)/) ? (args =~ /-k\s+(\d+)/)[0][1] : '24'
    // -R: memory for the reads, in MiB (Nonpareil's default 1024 caps it otherwise); 90 % of the task's
    def ram_mib = task.memory ? Math.min((task.memory.toMega() * 0.9) as long, 4194303L) : 1024
    """
    # R1 reads of at least k bp (nonpareil stops on a shorter one); one record per line, so the
    # image's busybox awk tests whole records (3.6x faster than four lines at a time: 52 s for 10 M reads)
    gzip -cdf ${r1} \\
        | paste - - - - \\
        | awk -F '\\t' -v k=${kmer} 'length(\$2) >= k { print; next } { short++ } END { print short + 0 > "reads_shorter_than_k.txt" }' \\
        | tr '\\t' '\\n' > ${prefix}_R1.fastq

    nonpareil \\
        -s ${prefix}_R1.fastq \\
        -T kmer \\
        -f fastq \\
        -b ${prefix} \\
        -t ${task.cpus} \\
        -R ${ram_mib} \\
        ${args}

    rm ${prefix}_R1.fastq

    nonpareil_fit.R ${prefix}.npo ${prefix} \$(cat reads_shorter_than_k.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(nonpareil -V 2>&1 | sed -n 's/^Nonpareil v//p')
        r-nonpareil: \$(Rscript -e 'cat(as.character(packageVersion("Nonpareil")))')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    printf 'sample\\treads\\treads_shorter_than_k\\tread_length\\tal\\tkappa\\tcoverage\\tdiversity\\tlr_star\\tdepth_95pct_pairs\\tmodel_r\\thas_model\\tconsistent\\twarning\\n' > ${prefix}_nonpareil.tsv
    printf 'sample\\tkind\\tx_adj\\tx_pairs\\tcoverage\\tcoverage_p25\\tcoverage_p75\\n' > ${prefix}_nonpareil_curve.tsv
    touch ${prefix}.npo ${prefix}.npl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: 3.5.5
        r-nonpareil: 3.5.3
    END_VERSIONS
    """
}
