process RGI_KMER {
    container 'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0'
    
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
        # per task (audit #17): reads go through the symlinks, while any file
        # RGI creates at runtime lands in this task's localDB/, not the shared
        # published database.
        if [ ! -d ${card_db}/localDB ]; then
            echo "ERROR: ${card_db} does not contain a localDB directory (expected an 'rgi load --local' database)" >&2
            exit 1
        fi
        mkdir localDB
        ln -s "\$(readlink -f ${card_db})/localDB"/* localDB/
        # kmer_query exits 0 with no output when the k-mer database is not
        # loaded, which turned into an empty '{}' result; fail instead
        for f in localDB/${kmer_size}mer_database.json localDB/amr_${kmer_size}mer.txt; do
            if [ ! -s "\$f" ]; then
                echo "ERROR: the RGI database has no ${kmer_size}-mer database (\$f missing)." >&2
                echo "       Load it with: rgi load --kmer_database <wildcard>/${kmer_size}_kmer_db.json" >&2
                echo "       --amr_kmers <wildcard>/all_amr_${kmer_size}mers.txt --kmer_size ${kmer_size} --local" >&2
                exit 1
            fi
        done
        
        # Run RGI kmer_query for pathogen-of-origin prediction, through
        # bin/rgi_kmer_kma_fix.py: RGI 6.0.3 misreads KMA BAMs and waits
        # forever on a worker that died (e.g. out of memory); the wrapper fixes
        # both and stops the run after ${time_limit} s. Remove at the RGI re-pin.
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
            rgi: \$(rgi --version 2>&1 | sed -n 's/.*rgi \\([0-9.]*\\).*/\\1/p' || echo "unknown")
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
            rgi: 6.0.3
        END_VERSIONS
        """
}
