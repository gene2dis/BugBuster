process VAMB {
    container 'quay.io/biocontainers/vamb:3.0.2--py36h91eb985_2'

    label 'process_medium'

    input:
        tuple val(meta), path(contigs), path(depth)
    
    output:
        tuple val(meta), path("*_vamb_bins"), emit: vamb

    script:
        def prefix = "${meta.id}"
        def min_fasta = "${params.vamb_min_fasta}"

        """
	vamb bin default \\
             --fasta ${contigs} \\
             --jgi ${depth} \\
             --outdir ${prefix}_vamb_bins \\
             --minfasta ${min_fasta}
	"""
}

process VAMB_COASSEMBLY {
    container 'quay.io/biocontainers/vamb:3.0.2--py36h91eb985_2'

    label 'process_high'

    input:
        path(bams_and_contigs)
    
    output:
        path("co_assembly_semibin_bins"), emit: semibin

    script:
        """
        contigs=`ls | grep -E '.+?filtered_contigs.+' | tr '\\n' ',' | sed 's/.\$//'`
        bam_list=`ls | grep -E '.+?all_reads.bam' | tr '\\n' ' ' | sed 's/.\$//'`


	"""
}
