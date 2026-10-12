process SYLPH_PROFILE {
    tag "${meta.id}"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/sylph:1.0.0--hb42e459_0' :
        'quay.io/biocontainers/sylph:1.0.0--hb42e459_0' }"

    label 'process_medium'

    input:
        tuple val(meta), path(reads)
        path(db)
        val(db_name)

    output:
        tuple val(meta), path("*_sylph_${db_name}.tsv"), emit: profile
        path("*_sylph_${db_name}_report.tsv")         , emit: report
        path "versions.yml"                            , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        // Paired reads only: sylph treats every -r file as a separate sample,
        // so a staged singleton cannot join the pair's profile (as for kraken2)
        """
        # -u scales Sequence_abundance by the estimated unknown fraction, so it
        # becomes each genome's share of ALL input reads: its sum is the
        # classified fraction and, times the read pairs, the estimated pairs
        # per genome (TAXONOMY_PHYLOSEQ). The header-only TSV sylph writes for
        # a sample without hits is handled downstream.
        sylph profile \\
            ${db} \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -u \\
            -t ${task.cpus} \\
            ${args} \\
            -o ${prefix}_sylph_${db_name}.tsv

        # Uncompressed .fastq/.fq reach this step with --quality_control
        # false; the image's busybox zcat has no -f passthrough, so pick the
        # reader from the gzip magic. A line count that is not a multiple of 4
        # is a truncated or malformed FASTQ, not a read count
        if [ "\$(head -c 2 ${reads[0]} | od -An -tx1 | tr -d ' \\n')" = "1f8b" ]; then
            reader=zcat
        else
            reader=cat
        fi
        read_pairs=\$(\${reader} ${reads[0]} | awk '
            END { if (NR % 4) { print "ERROR: " NR " lines in ${reads[0]}, not a multiple of 4 (truncated or malformed FASTQ)" > "/dev/stderr"; exit 1 }
                  printf "%.0f", NR / 4 }')

        classified=\$(awk -F'\\t' '
            NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "Sequence_abundance") col = i
                      if (!col) { print "ERROR: no Sequence_abundance column in sylph output" > "/dev/stderr"; exit 1 }
                      next }
            { sum += \$col }
            END { v = sum / 100; if (v > 1) v = 1; printf "%.6f", v }' ${prefix}_sylph_${db_name}.tsv)
        unclassified=\$(awk -v c="\${classified}" 'BEGIN { printf "%.6f", 1 - c }')

        printf 'Id\\tSylph DB\\tUnclassified\\tClassified\\tRead pairs\\n' > ${prefix}_sylph_${db_name}_report.tsv
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${prefix}" "${db_name}" "\${unclassified}" "\${classified}" "\${read_pairs}" \\
            >> ${prefix}_sylph_${db_name}_report.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sylph: \$(sylph --version 2>&1 | sed 's/sylph //')
            sylph_db: ${db.name}
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        printf 'Sample_file\\tGenome_file\\tTaxonomic_abundance\\tSequence_abundance\\tAdjusted_ANI\\tTrue_cov\\tANI_5-95_percentile\\tEff_lambda\\tLambda_5-95_percentile\\tMedian_cov\\tMean_cov_geq1\\tContainment_ind\\tNaive_ANI\\tkmers_reassigned\\tContig_name\\n' \\
            > ${prefix}_sylph_${db_name}.tsv
        printf 'Id\\tSylph DB\\tUnclassified\\tClassified\\tRead pairs\\n' > ${prefix}_sylph_${db_name}_report.tsv
        printf '%s\\t%s\\t1.000000\\t0.000000\\t0\\n' "${prefix}" "${db_name}" >> ${prefix}_sylph_${db_name}_report.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sylph: 1.0.0
            sylph_db: ${db.name}
        END_VERSIONS
        """
}
