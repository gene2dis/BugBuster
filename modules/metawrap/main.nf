process METAWRAP {
    container 'quay.io/biocontainers/metawrap:1.2--hdfd78af_2'

    cpus 10

    publishDir "${params.output}/workflow/${meta.id}/metawrap", pattern: '*metawrap*bins'
    publishDir "${params.output}/Assembly/${meta.id}/Refined_bins", mode: 'copy', pattern: '*metawrap*bins'

    input:
        tuple val(meta), path(metabat2_bins), path(semibin_bins), path(autometa_bins)
	path(metawrap_db)

    output:
        tuple val(meta), path("*metawrap*bins"), emit: metawrap

    cache 'lenient'

    script:
        def prefix = "${meta.id}"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
	echo ${metawrap_db} | checkm data setRoot ${metawrap_db}
        cp -rL ${metabat2_bins} metabat_wp_bins
	cp -rL ${semibin_bins} semibin_wp_bins
        cp -rL ${autometa_bins} autometa_wp_bins
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
