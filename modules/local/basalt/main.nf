process METAWRAP {
    container 'quay.io/ffuentessantander/basalt:v1.1.0'

    label 'process_high'

    input:
        tuple val(meta), path(metabat2_bins), path(semibin_bins), path(comebin_bins)
	path(metawrap_db)

    output:
        tuple val(meta), path("*metawrap*bins"), emit: metawrap

    script:
        def prefix = "${meta.id}"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
	echo ${metawrap_db} | checkm data setRoot ${metawrap_db}
        cp -r ${metabat2_bins} metabat_wp_bins
	cp -r ${semibin_bins} semibin_wp_bins
        cp -r ${comebin_bins} comebin_wp_bins
	metawrap bin_refinement \\
		-o Refined_bins \\
		-t $task.cpus \\
		-A metabat_wp_bins \\
		-B semibin_wp_bins \\
		-C comebin_wp_bins \\
		-c ${completeness} \\
                -x ${contamination}
	chmod 777 Refined_bins
	mv Refined_bins/metawrap_${completeness}_${contamination}_bins ./${prefix}_metawrap_${completeness}_${contamination}_bins
        rm -r metabat_wp_bins
        rm -r semibin_wp_bins
        rm -r comebin_wp_bins
	"""
}

process METAWRAP_COASSEMBLY {
    container 'quay.io/biocontainers/metawrap:1.2--hdfd78af_2'

    label 'process_high'

    input:
        path(bins)
	path(metawrap_db)

    output:
        path("*metawrap*bins"), emit: metawrap

    script:
        def prefix = "co_assembly"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
	echo ${metawrap_db} | checkm data setRoot ${metawrap_db}
        cp -r ${bins[0]} metabat_wp_bins
	cp -r ${bins[1]} semibin_wp_bins
        cp -r ${bins[2]} autometa_wp_bins
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
