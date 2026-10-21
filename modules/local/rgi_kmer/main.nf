process RGI_KMER {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.8--pyh05cac1d_0' :
        'quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0' }"
    
    label 'process_medium'
    
    tag "${meta.id}"

    input:
        tuple val(meta), path(bam)
        path(card_db)

    output:
        // rgi kmer_query --bwt writes <out>_<k>mer_analysis.{json,allele.txt,gene.txt}
        tuple val(meta), path("*_${params.rgi_kmer_size ?: 61}mer_analysis.json"), emit: kmer_json
        tuple val(meta), path("*.allele.txt"), emit: allele_predictions, optional: true
        tuple val(meta), path("*.gene.txt"), emit: gene_predictions, optional: true
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        def kmer_size = params.rgi_kmer_size ?: 61
        def minimum = params.rgi_min_kmer_coverage ?: 10
        // Wall-clock cap enforced by the wrapper (modules.config ext.time_limit)
        def time_limit = ((task.ext.time_limit ?: 2.h) as Duration).toSeconds()
        
        """
        set +u
        if [ -f /usr/local/env-execute ]; then
            # env-execute resets PATH; keep the pipeline's bin/ on it
            pipeline_path="\$PATH"
            source /usr/local/env-execute
            export PATH="\$PATH:\$pipeline_path"
        fi
        set -u

        # RGI --local looks for localDB/ in the current directory. Symlink the
        # DB contents per file instead of copying the multi-GB CARD/WildCARD DB
        # per task (audit #17): cp -rsL makes real directories (following
        # symlinked ones, localDB itself included) and a symlink per file, so
        # any file RGI creates at runtime lands in this task's localDB/, not
        # the shared published database (as in RGI_BWT; kmer_query only reads).
        if [ ! -d ${card_db}/localDB ]; then
            echo "ERROR: ${card_db} does not contain a localDB directory (expected an 'rgi load --local' database)" >&2
            exit 1
        fi
        cp -rsL "\$(readlink -f ${card_db}/localDB)" localDB
        # kmer_query exits 0 with no output when the k-mer database is not
        # loaded. Say so and skip it: the empty JSON keeps the declared output,
        # and the summary reports the sample's RGI hosts as 'no_kmer_data'
        kmer_db_ok=true
        for f in localDB/${kmer_size}mer_database.json localDB/amr_${kmer_size}mer.txt; do
            if [ ! -s "\$f" ]; then
                echo "WARNING: the RGI database has no ${kmer_size}-mer database (\$f missing);" >&2
                echo "         pathogen of origin is skipped for ${prefix}. Load it with: rgi load" >&2
                echo "         --kmer_database <wildcard>/${kmer_size}_kmer_db.json --amr_kmers" >&2
                echo "         <wildcard>/all_amr_${kmer_size}mers.txt --kmer_size ${kmer_size} --local" >&2
                kmer_db_ok=false
            fi
        done
        if [ "\$kmer_db_ok" = false ]; then
            # a new file, never written through a staged input of the same name
            rm -f ${prefix}_rgi_kmer_${kmer_size}mer_analysis.json
            echo '{}' > ${prefix}_rgi_kmer_${kmer_size}mer_analysis.json
            printf '"%s":\\n    rgi: %s\\n' "${task.process}" "\$(rgi main --version 2>/dev/null || echo unknown)" > versions.yml
            exit 0
        fi
        
        # Run RGI kmer_query for pathogen-of-origin prediction, through
        # bin/rgi_kmer_kma_fix.py: RGI waits forever on a worker that died
        # (e.g. out of memory); the wrapper fails the task instead and stops
        # the run after ${time_limit} s. Remove it when RGI checks its workers.
        # Every worker is a fork of the process holding the k-mer database, so
        # memory grows with -n (task.cpus, set low in modules.config).
        rgi_kmer_kma_fix.py --time-limit ${time_limit} kmer_query \\
            --bwt \\
            -i ${bam} \\
            -k ${kmer_size} \\
            -n ${task.cpus} \\
            -m ${minimum} \\
            -o ${prefix}_rgi_kmer \\
            --local \\
            ${args}

        # Verify output files were created
        if [ ! -f "${prefix}_rgi_kmer_${kmer_size}mer_analysis.json" ]; then
            echo "WARNING: RGI kmer_query did not produce expected JSON output"
            # Create empty output to prevent pipeline failure
            echo '{}' > ${prefix}_rgi_kmer_${kmer_size}mer_analysis.json
        fi

        # Create versions file
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: \$(rgi main --version 2>/dev/null || echo "unknown")
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        def kmer_size = params.rgi_kmer_size ?: 61
        """
        touch ${prefix}_rgi_kmer_${kmer_size}mer_analysis.json
        touch ${prefix}_rgi_kmer_${kmer_size}mer_analysis.allele.txt
        touch ${prefix}_rgi_kmer_${kmer_size}mer_analysis.gene.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.8
        END_VERSIONS
        """
}
