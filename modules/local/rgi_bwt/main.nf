process RGI_BWT {
    container 'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0'
    
    label 'process_medium'
    
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)
        path(card_db)

    output:
        tuple val(meta), path("*.allele_mapping_data.txt"), emit: allele_mapping
        tuple val(meta), path("*.gene_mapping_data.txt"), emit: gene_mapping
        tuple val(meta), path("*.sorted.length_100.bam"), emit: bam, optional: true
        path("*.overall_mapping_stats.txt"), emit: overall_stats
        path("*.reference_mapping_stats.txt"), emit: reference_stats
        path("*.artifacts_mapping_stats.txt"), emit: artifacts_stats, optional: true
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        def aligner = params.rgi_aligner ?: 'kma'
        def wildcard = params.rgi_include_wildcard ? '--include_wildcard' : ''
        
        // rgi bwt only accepts --read_one/--read_two, so read-level RGI runs on
        // the paired reads only; singleton reads (reads[2], carried since the
        // audit #9 QC rework) are deliberately excluded (documented in the manual)
        def read_one = reads[0]
        def read_two = reads.size() > 1 ? reads[1] : ''
        
        """
        set +u
        if [ -f /usr/local/env-execute ]; then
            source /usr/local/env-execute
        fi
        set -u

        # Copy database to local directory for RGI access
        cp -r ${card_db} rgi_db
        
        # RGI --local flag looks for localDB/ in current directory
        ln -s rgi_db/localDB localDB
        
        # Run RGI bwt (stderr redirected to log file to avoid flooding with low-coverage warnings)
        rgi bwt \\
            --read_one ${read_one} \\
            ${read_two ? "--read_two ${read_two}" : ""} \\
            --aligner ${aligner} \\
            --threads ${task.cpus} \\
            --output_file ${prefix}_rgi_bwt \\
            --local \\
            ${wildcard} \\
            ${args} \\
            2>${prefix}_rgi_bwt.log

        # Low-coverage fallback (audit #8): if RGI produced no mapping results,
        # create every declared non-optional output (NA allele/gene tables plus
        # both stats files). No BAM is created - the bam emit is optional and
        # RGI_KMER is skipped for this sample (no alignments to query).
        if [ ! -f "${prefix}_rgi_bwt.allele_mapping_data.txt" ]; then
            echo "WARNING: RGI bwt did not produce expected output files - creating dummy outputs for low coverage sample"
            rgi_bwt_lowcov_fallback.sh ${prefix}_rgi_bwt
        fi

        # Create versions file
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: \$(rgi --version 2>&1 | sed -n 's/.*rgi \\([0-9.]*\\).*/\\1/p' || echo "unknown")
            ${aligner}: \$(${aligner} --version 2>&1 | head -1 || echo "unknown")
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        touch ${prefix}_rgi_bwt.allele_mapping_data.txt
        touch ${prefix}_rgi_bwt.gene_mapping_data.txt
        touch ${prefix}_rgi_bwt.sorted.length_100.bam
        touch ${prefix}_rgi_bwt.overall_mapping_stats.txt
        touch ${prefix}_rgi_bwt.reference_mapping_stats.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.3
            kma: 1.4.14
        END_VERSIONS
        """
}
