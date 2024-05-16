process CALCULATE_DEPTH {
    container 'quay.io/biocontainers/metabat2:2.15--h4da6f23_2'

    cpus 5

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
