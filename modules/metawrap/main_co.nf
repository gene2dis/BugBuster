process METAWRAP {
    container 'quay.io/biocontainers/metawrap:1.2--hdfd78af_2'

    cpus 40

    publishDir "${params.output}/workflow/co_assembly/metawrap", pattern: '*metawrap*bins'
    publishDir "${params.output}/Co_assembly/Refined_bins", mode: 'copy', pattern: '*metawrap*bins'

    input:
        path(bins)
	path(metawrap_db)

    output:
        path("*metawrap*bins"), emit: metawrap

    cache 'lenient'

    script:
        def prefix = "co_assembly"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
	echo ${metawrap_db} | checkm data setRoot ${metawrap_db}
        cp -rL ${bins[0]} metabat_wp_bins
	cp -rL ${bins[1]} semibin_wp_bins
        cp -rL ${bins[2]} autometa_wp_bins
	metawrap bin_refinement \\
		-o Refined_bins \\
		-t $task.cpus \\
		-A metabat_wp_bins \\
		-B semibin_wp_bins \\
		-C autometa_wp_bins \\
		-c ${completeness} \\
                -x ${contamination}
	chmod 777 Refined_bins
	mv Refined_bins/metawrap_${completeness}_${contamination}_bins ./${prefix}_metawrap_${completeness}_${contamination}_bins
	"""
}
