process COMEBIN {
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_0'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/comebin", pattern: '*_comebin_bins/comebin_res/comebin_res_bins'
    publishDir "${params.output}/Assembly/${meta.id}/Bins", mode: 'copy', pattern: '*_comebin_bins/comebin_res/comebin_res_bins'
    publishDir "${params.output}/Assembly/${meta.id}/comebin_emmbeding", mode: 'copy', pattern: '*_comebin_bins/comebin_res/embeddings.tsv'
    publishDir "${params.output}/Assembly/${meta.id}/comebin_emmbeding", mode: 'copy', pattern: '*_comebin_bins/comebin_res/covembeddings.tsv'

    input:
        tuple val(meta), path(contigs), path(bam)
    
    output:
        tuple val(meta), path("*_comebin_bins/comebin_res/comebin_res_bins"), emit: comebin

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
	"""
}
