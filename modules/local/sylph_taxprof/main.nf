process SYLPH_TAXPROF {
    tag "${meta.id}"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/sylph-tax:1.9.2--pyhdfd78af_0' :
        'quay.io/biocontainers/sylph-tax:1.9.2--pyhdfd78af_0' }"

    label 'process_single'

    input:
        tuple val(meta), path(profile)
        path(taxonomy)
        val(db_name)

    output:
        tuple val(meta), path("${meta.id}_sylph_${db_name}.sylphmpa"), emit: sylphmpa
        path "versions.yml"                                          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        def header = 'clade_name\\trelative_abundance\\tsequence_abundance\\tANI (if strain-level)\\tCoverage (if strain-level)'
        // --no-config (a global option, before the subcommand) with the staged
        // metadata file as -t: no ~/.config/sylph-tax and no downloaded
        // taxonomy directory, so the read-only images and HOME-less HPC
        // nodes work
        """
        if [ "\$(wc -l < ${profile})" -le 1 ]; then
            # No genome passed sylph's threshold: sylph-tax writes no file for
            # a header-only profile, so write its two header lines ourselves
            printf '#SampleID\\t%s\\tTaxonomies_used:[%s]\\n${header}\\n' "${prefix}" "'${taxonomy}'" \\
                > ${prefix}_sylph_${db_name}.sylphmpa
        else
            sylph-tax --no-config taxprof \\
                ${profile} \\
                -t ${taxonomy} \\
                ${args} \\
                -o ${prefix}_tmp_

            # One profile holds one sample, so exactly one .sylphmpa
            # (<prefix><basename of the sample's R1>.sylphmpa)
            set -- ${prefix}_tmp_*.sylphmpa
            if [ "\$#" -ne 1 ] || [ ! -f "\$1" ]; then
                echo "ERROR: expected one .sylphmpa from sylph-tax, found \$#" >&2
                exit 1
            fi
            mv "\$1" ${prefix}_sylph_${db_name}.sylphmpa
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sylph-tax: \$(sylph-tax --version 2>&1 | sed 's/sylph-tax //')
            sylph_taxonomy: ${taxonomy.name}
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        printf '#SampleID\\t%s\\tTaxonomies_used:[stub]\\nclade_name\\trelative_abundance\\tsequence_abundance\\tANI (if strain-level)\\tCoverage (if strain-level)\\n' "${prefix}" \\
            > ${prefix}_sylph_${db_name}.sylphmpa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sylph-tax: 1.9.2
            sylph_taxonomy: ${taxonomy.name}
        END_VERSIONS
        """
}
