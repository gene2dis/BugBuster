process PUBLISH_CHECKM2_REPORTS {
    tag "meta.id"
    label 'process_single'
    
    publishDir path: { "${params.output}/04_binning/per_sample/${meta.id}/quality/checkm2" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reports)

    output:
    path "*.tsv", emit: reports

    script:
    """
    # Files are already staged by Nextflow, publishDir handles publishing
    # List files for debugging
    ls -la *.tsv 2>/dev/null || echo "No TSV files found"
    """
}
