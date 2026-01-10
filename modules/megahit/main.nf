/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MEGAHIT Assembly Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenome assembly using MEGAHIT
----------------------------------------------------------------------------------------
*/

process MEGAHIT {
    tag "${meta.id}"
    container 'quay.io/biocontainers/megahit:1.2.9--h43eeafb_4'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/assembly", mode: params.publish_dir_mode, pattern: "*_contigs.fa"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), path("*_contigs.fa"), emit: contigs_and_reads
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def r1 = reads[0]
    def r2 = reads[1]
    def singleton = reads.size() > 2 ? reads[2] : null
    def singleton_arg = singleton ? "--read ${singleton}" : ""

    """
    megahit \\
        -1 ${r1} \\
        -2 ${r2} \\
        ${singleton_arg} \\
        -t ${task.cpus} \\
        -o ${prefix}_assembly

    mv ${prefix}_assembly/final.contigs.fa ${prefix}_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(megahit --version 2>&1 | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: 1.2.9
    END_VERSIONS
    """
}

process MEGAHIT_COASSEMBLY {
    tag "coassembly"
    container 'quay.io/biocontainers/megahit:1.2.9--h43eeafb_4'

    label 'process_high'

    publishDir "${params.output}/workflow/coassembly/assembly", mode: params.publish_dir_mode, pattern: "*_contigs.fa"

    input:
    path(reads)

    output:
    tuple path(reads), path("co_assembly_contigs.fa"), emit: megahit
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Build comma-separated file lists for MEGAHIT
    R1_list=\$(ls -1 | grep -E '_R1_|_1\\.' | tr '\\n' ',' | sed 's/,\$//')
    R2_list=\$(ls -1 | grep -E '_R2_|_2\\.' | tr '\\n' ',' | sed 's/,\$//')
    Singleton_list=\$(ls -1 | grep -iE 'singleton|unpaired' | tr '\\n' ',' | sed 's/,\$//')

    if [ -n "\$Singleton_list" ]; then
        megahit \\
            -1 \$R1_list \\
            -2 \$R2_list \\
            --read \$Singleton_list \\
            -t ${task.cpus} \\
            -o co_assembly
    else
        megahit \\
            -1 \$R1_list \\
            -2 \$R2_list \\
            -t ${task.cpus} \\
            -o co_assembly
    fi

    mv co_assembly/final.contigs.fa co_assembly_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(megahit --version 2>&1 | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    """
    touch co_assembly_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: 1.2.9
    END_VERSIONS
    """
}
