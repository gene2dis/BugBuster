// Clusters predicted ARG proteins at four fixed identity tiers (90/95/99/100%).
// The tiers are deliberate: the former single-valued mmseqs_* params could not
// describe this design and were removed (audit #14).
process CLUSTERING {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0' :
        'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0' }"

    tag "arg_clustering"
    label 'process_high'

    input:
        path(arg_faa)

    output:
        path("*_cluster.tsv"), emit: clusters
        path "versions.yml", emit: versions

    script:

        """
        cat *.faa > mmseq_db.faa
        mmseqs easy-cluster mmseq_db.faa Arg_cluster_90 tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0 --threads $task.cpus
        mmseqs easy-cluster mmseq_db.faa Arg_cluster_95 tmp --min-seq-id 0.95 -c 0.8 --cov-mode 0 --threads $task.cpus
        mmseqs easy-cluster mmseq_db.faa Arg_cluster_99 tmp --min-seq-id 0.99 -c 0.9 --cov-mode 0 --threads $task.cpus
        mmseqs easy-cluster mmseq_db.faa Arg_cluster_100 tmp --min-seq-id 1.0 -c 0.9 --cov-mode 0 --threads $task.cpus
        rm mmseq_db.faa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mmseqs: \$(mmseqs version)
        END_VERSIONS
	"""

    stub:
        """
        touch Arg_cluster_90_cluster.tsv
        touch Arg_cluster_95_cluster.tsv
        touch Arg_cluster_99_cluster.tsv
        touch Arg_cluster_100_cluster.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mmseqs: 18.8cc5c
        END_VERSIONS
        """
}
