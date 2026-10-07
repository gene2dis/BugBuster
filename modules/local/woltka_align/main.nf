/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WOLTKA_ALIGN Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Bowtie2 alignment of host-removed reads against the Web of Life (WoLr2)
    genome database, the first half of the Woltka read-level functional
    backend (design doc Section 4.6.2, task T8a). Split from WOLTKA_CLASSIFY
    so the expensive, memory-heavy alignment (the WoLr2 index needs >= 68 GB
    RAM) retries and is resourced independently of the cheap classification
    (eggNOG search/annotate precedent).

    Alignment flags come from task.ext.args (config/modules.config): the
    SHOGUN multi-hit set from the WoLr2 Bowtie2 README by default (owner
    decision 2026-10-07; Woltka then divides each read 1/k across its k
    hits, unless --woltka_uniq). --no-head --no-unal are fixed here because
    the output contract depends on them. Reads are aligned PAIRED (-1/-2)
    plus -U for the singleton file: Woltka counts mates as separate queries
    by SAM flag, which only works when the mates carry pair flags.

    The SAM is trimmed as the WoLr2 README recommends (columns 1-9, SEQ and
    QUAL replaced by '*') and gzipped - Woltka only needs positions, and
    this is what keeps the multi-hit intermediate small (Section 4.6.2's
    scratch-space note). The index basename is discovered under
    <db>/databases/bowtie2/ (WoLr2 ships .bt2l large-index files; .bt2 is
    accepted too).

    Input:  tuple val(meta), path(reads), path(wol_db)
            reads = [R1, R2] or [R1, R2, Singleton] (or a single file)
    Output: tuple val(meta), path("*.wol.sam.gz"), emit: sam
            tuple val(meta), path("*.bowtie2.log"), emit: log
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process WOLTKA_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.3--py310ha0a81b8_0' :
        'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0' }"

    input:
    tuple val(meta), path(reads), path(wol_db)

    output:
    tuple val(meta), path("*.wol.sam.gz") , emit: sam
    tuple val(meta), path("*.bowtie2.log"), emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_list = reads instanceof List ? reads : [reads]
    def read_args = read_list.size() >= 2
        ? "-1 ${read_list[0]} -2 ${read_list[1]}" + (read_list.size() > 2 ? " -U ${read_list[2]}" : '')
        : "-U ${read_list[0]}"
    """
    index_file=\$(ls ${wol_db}/databases/bowtie2/*.rev.1.bt2l ${wol_db}/databases/bowtie2/*.rev.1.bt2 2>/dev/null | head -1 || true)
    if [ -z "\${index_file}" ]; then
        echo "ERROR: no Bowtie2 index (*.rev.1.bt2l / *.rev.1.bt2) under ${wol_db}/databases/bowtie2/ - expected the WoLr2 layout" >&2
        exit 1
    fi
    index="\${index_file%.rev.1.bt2*}"

    bowtie2 \\
        -p ${task.cpus} \\
        -x "\${index}" \\
        ${read_args} \\
        ${args} \\
        --no-head \\
        --no-unal \\
        2> ${prefix}.bowtie2.log \\
        | cut -f1-9 \\
        | sed 's/\$/\\t*\\t*/' \\
        | gzip > ${prefix}.wol.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "" | gzip > ${prefix}.wol.sam.gz
    touch ${prefix}.bowtie2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.3
    END_VERSIONS
    """
}
