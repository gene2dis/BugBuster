process CALCULATE_DEPTH_COASSEMBLY {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    cpus 1

    publishDir "${params.output}/workflow/co_assembly/depth_file", pattern: '*depth.txt'

    input:
       path(bam)

    output:
       path("*depth.txt"), emit: depth

    script:

        """
        bam_list=`ls | grep -E '.+?all_reads.bam' | tr '\\n' ' ' | sed 's/.\$//'`
	jgi_summarize_bam_contig_depths --outputDepth depth.txt \${bam_list}
	mv depth.txt co_assembly_depth.txt
	"""
}

process CALCULATE_DEPTH {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    cpus 1

    publishDir "${params.output}/workflow/${meta.id}/depth_file", pattern: '*depth.txt'

    input:
       tuple val(meta), path(contigs), path(bam)

    output:
       tuple val(meta), path(contigs), path("*depth.txt"), emit: depth

    script:
        def prefix = "${meta.id}"

        """
	jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
	mv depth.txt ${prefix}_depth.txt
	"""
}
