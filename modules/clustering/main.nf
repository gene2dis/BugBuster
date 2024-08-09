process CLUSTERING {
    container 'quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_1'

    label 'process_medium'

    input:
        path(mmseq_results)

    output:
        path("*_cluster_*"), emit: clusters

    script:

        """
        mmseqs easy-cluster ${mmseq_results} Arg_cluster_90 tmp --min-seq-id 0.9 -c 0.8 --cov-mode 0 --threads $task.cpus
        mmseqs easy-cluster ${mmseq_results} Arg_cluster_95 tmp --min-seq-id 0.95 -c 0.8 --cov-mode 0 --threads $task.cpus
        mmseqs easy-cluster ${mmseq_results} Arg_cluster_99 tmp --min-seq-id 0.99 -c 0.9 --cov-mode 0 --threads $task.cpus
        mmseqs easy-cluster ${mmseq_results} Arg_cluster_100 tmp --min-seq-id 1.0 -c 0.9 --cov-mode 0 --threads $task.cpus
	"""
}
