/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BBMAP Contig Filtering Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Filter assembled contigs by minimum length using BBMap reformat.sh
    
    Supports both per-sample and co-assembly modes
    
    Input:
        tuple val(meta), path(reads), path(contigs)
    
    Output:
        contigs_with_reads: tuple val(meta), path(reads), path(filtered_contigs)
        contigs_only: tuple val(meta), path(filtered_contigs)
        stats: path(contig.stats)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process BBMAP {
    tag "${meta.id}"
    label 'process_single'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.06--h92535d8_0' :
        'quay.io/biocontainers/bbmap:39.06--h92535d8_0' }"

    // Use storeDir for filtered contigs to enable immediate work dir cleanup
    // This stores outputs permanently and cleans work directory automatically
    storeDir params.store_filtered_contigs ? "${params.output}/assembly/${meta.id}" : null

    input:
    tuple val(meta), path(reads), path(contigs)

    output:
    tuple val(meta), path(reads), path("${meta.id}_filtered_contigs.fa"), emit: contigs_with_reads
    tuple val(meta), path("${meta.id}_filtered_contigs.fa")            , emit: contigs_only
    path "${meta.id}_contig.stats"                                     , emit: stats
    path "${meta.id}_filter_report.txt"                                , emit: filter_report
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    def min_length = task.ext.min_length ?: params.bbmap_length ?: 1000
    
    // Input validation
    if (!meta || !meta.id) {
        error "BBMAP: meta.id is required"
    }
    if (!contigs) {
        error "BBMAP: No contig file provided for sample ${meta.id}"
    }
    if (min_length < 0) {
        error "BBMAP: Minimum length must be non-negative (got: ${min_length})"
    }

    """
    set -euo pipefail
    
    # Validate contig file exists
    if [ ! -f "${contigs}" ]; then
        echo "ERROR: Contig file not found: ${contigs}" >&2
        exit 1
    fi
    
    # Count contigs before filtering
    contigs_before=\$(grep -c '^>' "${contigs}" || echo 0)
    
    # BBMap reformat.sh handles both gzipped and uncompressed contigs
    reformat.sh \\
        in=${contigs} \\
        out=${prefix}_filtered_contigs.fa \\
        minlength=${min_length} \\
        ${args}
    
    # Count contigs after filtering
    contigs_after=\$(grep -c '^>' ${prefix}_filtered_contigs.fa || echo 0)
    
    # Generate filter report
    cat > ${prefix}_filter_report.txt <<-REPORT
sample_id\t${prefix}
contigs_before_filter\t\${contigs_before}
contigs_after_filter\t\${contigs_after}
min_length_threshold\t${min_length}
all_contigs_filtered\t\$([ \${contigs_after} -eq 0 ] && echo "TRUE" || echo "FALSE")
REPORT
    
    # Log warning if all contigs were filtered
    if [ \${contigs_after} -eq 0 ]; then
        echo "WARNING: All \${contigs_before} contigs filtered out for sample ${prefix} (min_length=${min_length})" >&2
    fi
    
    stats.sh \\
        in=${contigs} \\
        out=${prefix}_contig.stats
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh 2>&1 | grep -o 'BBMap version.*' | sed 's/BBMap version //')
    END_VERSIONS
    """

    stub:
    prefix = meta.id
    """
    # Create a minimal filtered contig file
    cat > ${prefix}_filtered_contigs.fa << 'EOF'
>filtered_contig_1 length=2000
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>filtered_contig_2 length=3500
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
EOF
    
    # Create a minimal stats file
    cat > ${prefix}_contig.stats << 'EOF'
#name	length	gc
 filtered_contig_1	2000	0.50
filtered_contig_2	3500	0.48
EOF
    
    # Create filter report
    cat > ${prefix}_filter_report.txt << 'EOF'
sample_id	${prefix}
contigs_before_filter	10
contigs_after_filter	2
min_length_threshold	1000
all_contigs_filtered	FALSE
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: 39.06
    END_VERSIONS
    """
}
