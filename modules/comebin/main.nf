process COMEBIN {
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_0'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/comebin", pattern: '*_comebin_res_bins'
    publishDir "${params.output}/Assembly/${meta.id}/Bins", mode: 'copy', pattern: '*_comebin_res_bins'
    publishDir "${params.output}/Assembly/${meta.id}/comebin_emmbeding", mode: 'copy', pattern: '*beddings.tsv'

    input:
        tuple val(meta), path(contigs), path(bam)
    
    output:
        tuple val(meta), path("*_comebin_res_bins"), emit: comebin

    script:
        def prefix = "${meta.id}"

        """
        mkdir ${prefix}_comebin_bins
        bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus

        mv ${prefix}_comebin_bins/comebin_res/comebin_res_bins ./${prefix}_comebin_res_bins
        mv ${prefix}_comebin_bins/comebin_res/embeddings.tsv ./${prefix}_comebin_embeddings.tsv
        mv ${prefix}_comebin_bins/comebin_res/covembeddings.tsv ./${prefix}_comebin_covembeddings.tsv

        rm -rf ${prefix}_comebin_bins
	"""
}
