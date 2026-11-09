/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTISAMPLE_CONCAT Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The assemblies of one --multisample_binning group concatenated into one
    reference, each contig named <sample>:<contig> (the first word of its
    header), each sequence on one line, the samples in the order given (sorted
    by id). Sample ids cannot contain ':' (INPUT_CHECK), so the prefix splits
    back unambiguously. The same file as SemiBin2 concatenate_fasta writes
    (roadmap Phase 18a, docs/benchmarks/multisample_binning_2026.md).

    Input:
        tuple val(meta), path(contigs)
            meta:    [id: group label, samples: [sample ids, sorted]]
            contigs: each sample's filtered contigs, in meta.samples order

    Output:
        fasta:    tuple val(meta), path("${meta.id}_multisample_contigs.fa")
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process MULTISAMPLE_CONCAT {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.18--h38e344b_2' :
        'quay.io/biocontainers/metabat2:2.18--h38e344b_2' }"

    input:
    tuple val(meta), path(contigs, stageAs: 'in?/*')

    output:
    tuple val(meta), path("${meta.id}_multisample_contigs.fa"), emit: fasta
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.id
    def files = contigs instanceof List ? contigs : [contigs]
    if (files.size() != meta.samples.size()) {
        error "MULTISAMPLE_CONCAT: ${files.size()} contig files for the ${meta.samples.size()} samples of group ${meta.id}"
    }
    def concat = [meta.samples, files].transpose().collect { sample, f -> "add_sample '${sample}' ${f}" }.join('\n    ')
    """
    set -euo pipefail

    # Header: <sample>:<first word>; the sequence on one line, printed as it
    # comes (a string grown line by line is quadratic in BusyBox awk)
    add_sample() {
        awk -v p="\$1:" '/^>/ { if (n++) printf "\\n"; name = substr(\$0, 2); sub(/[ \\t].*/, "", name); print ">" p name; next }
            { printf "%s", \$0 }
            END { if (n) printf "\\n" }' "\$2"
    }

    {
    ${concat}
    } > ${prefix}_multisample_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: \$(busybox 2>&1 | head -1 | sed -E 's/BusyBox v([0-9.]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    for s in ${meta.samples.join(' ')}; do
        printf '>%s:contig_1\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n' "\$s"
    done > ${prefix}_multisample_contigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: 1.36.1
    END_VERSIONS
    """
}
