/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTISAMPLE_SPLIT_DEPTH Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Splits the jgi_summarize_bam_contig_depths table of a --multisample_binning
    group's concatenation (one depth and variance column per sample of the
    group) into one table per assembly: the rows of the contigs named
    <sample>:<contig>, with the prefix removed, under the same header. MetaBAT2
    bins each assembly on its table. These depths are competitive (a read of a
    genome assembled in several samples counts on one copy only), so they are
    MetaBAT2's input only, never abundance.

    Input:
        tuple val(meta), path(depth)
            meta:  [id: group label, samples: [sample ids]]
            depth: CALCULATE_DEPTH's table of the concatenation

    Output:
        depth:    path("*_multisample_depth.txt") - <sample>_multisample_depth.txt
                  for every sample of the group (header only when its
                  assembly has no contigs)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process MULTISAMPLE_SPLIT_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.18--h38e344b_2' :
        'quay.io/biocontainers/metabat2:2.18--h38e344b_2' }"

    input:
    tuple val(meta), path(depth)

    output:
    path "*_multisample_depth.txt", emit: depth
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # Every sample gets a table with the header; each row goes to its
    # sample's table, the name without '<sample>:' (sample ids have no ':')
    awk -F '\\t' -v OFS='\\t' -v samples='${meta.samples.join(' ')}' '
        NR == 1 {
            n = split(samples, s, " ")
            for (i = 1; i <= n; i++) { want[s[i]] = 1; print > (s[i] "_multisample_depth.txt") }
            next
        }
        {
            c = index(\$1, ":")
            id = substr(\$1, 1, c - 1)
            if (c == 0 || !(id in want)) { print "ERROR: contig " \$1 " is named after no sample of the group" > "/dev/stderr"; exit 1 }
            \$1 = substr(\$1, c + 1)
            print > (id "_multisample_depth.txt")
        }' ${depth}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: \$(busybox 2>&1 | head -1 | sed -E 's/BusyBox v([0-9.]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    """
    for s in ${meta.samples.join(' ')}; do
        head -n 1 ${depth} > "\${s}_multisample_depth.txt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: 1.36.1
    END_VERSIONS
    """
}
