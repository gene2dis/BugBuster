process METAWRAP {
    container 'quay.io/ffuentessantander/metawrap:1.2'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/metawrap", pattern: '*metawrap*bins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/refined_bins", mode: 'copy', pattern: '*metawrap*bins'

    input:
        tuple val(meta), path(metabat2_bins), path(semibin_bins), path(comebin_bins)

    output:
        tuple val(meta), path("*metawrap*bins"), emit: metawrap

    cache 'lenient'

    script:
        def prefix = "${meta.id}"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
        cp -rL ${metabat2_bins} metabat_wp_bins
	cp -rL ${semibin_bins} semibin_wp_bins
        cp -rL ${comebin_bins} comebin_wp_bins
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
    container 'quay.io/ffuentessantander/metawrap:1.2'

    label 'process_high'

    publishDir "${params.output}/workflow/co_assembly/metawrap", pattern: '*metawrap*bins'
    publishDir "${params.output}/contigs_and_bins/co_assembly/refined_bins", mode: 'copy', pattern: '*metawrap*bins'

    input:
        path(bins)

    output:
        path("*metawrap*bins"), emit: metawrap

    cache 'lenient'

    script:
        def prefix = "co_assembly"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
        cp -rL ${bins[0]} metabat_wp_bins
	cp -rL ${bins[1]} semibin_wp_bins
        cp -rL ${bins[2]} comebin_wp_bins
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
