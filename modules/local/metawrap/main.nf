process METAWRAP {
    container 'quay.io/ffuentessantander/metawrap:1.2'

    label 'process_high'

    input:
        tuple val(meta), path(metabat2_bins), path(semibin_bins), path(comebin_bins)

    output:
        tuple val(meta), path("*metawrap*bins"), emit: metawrap

    script:
        def prefix = "${meta.id}"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
        mkdir -p metabat_wp_bins semibin_wp_bins comebin_wp_bins
        
        for bin in ${metabat2_bins}; do
            if [[ "\$bin" == *.fa.gz ]]; then
                gunzip -c "\$bin" > "metabat_wp_bins/\$(basename "\$bin" .gz)"
            elif [[ "\$bin" == *.fa ]]; then
                cp "\$bin" metabat_wp_bins/
            fi
        done
        
        cp -rL ${semibin_bins}/* semibin_wp_bins/ 2>/dev/null || true
        cp -rL ${comebin_bins}/* comebin_wp_bins/ 2>/dev/null || true
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

    input:
        path(bins)

    output:
        path("*metawrap*bins"), emit: metawrap

    script:
        def prefix = "co_assembly"
        def completeness = "${params.metawrap_completeness}"
        def contamination = "${params.metawrap_contamination}"

        """
        mkdir -p metabat_wp_bins semibin_wp_bins comebin_wp_bins
        
        for bin in ${bins[0]}; do
            if [[ "\$bin" == *.fa.gz ]]; then
                gunzip -c "\$bin" > "metabat_wp_bins/\$(basename "\$bin" .gz)"
            elif [[ "\$bin" == *.fa ]]; then
                cp "\$bin" metabat_wp_bins/
            fi
        done
        
        cp -rL ${bins[1]}/* semibin_wp_bins/ 2>/dev/null || true
        cp -rL ${bins[2]}/* comebin_wp_bins/ 2>/dev/null || true
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
