process COMEBIN {
    container 'quay.io/biocontainers/comebin:1.0.4--hdfd78af_0'

    label 'process_high'

    publishDir "${params.output}/workflow/${meta.id}/comebin", pattern: '*_comebin_bins/comebin_res/comebin_res_bins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/raw_bins", mode: 'copy', pattern: '*_comebin_bins/comebin_res/comebin_res_bins'

    input:
        tuple val(meta), path(contigs), path(bam)
    
    output:
        tuple val(meta), path("*_comebin_bins/comebin_res/comebin_res_bins"), emit: comebin

    script:
        def prefix = "${meta.id}"

        """
        mkdir ${prefix}_comebin_bins

        if [[ $task.attempt == 1 ]]; then
               bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus 
        elif [[ $task.attempt == 2 ]]; then 
                bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus \\
                   -b 896 \\
                   -e 1792 \\
                   -c 1792
        elif [[ $task.attempt == 3 ]]; then
                bash run_comebin.sh \\
                   -a ${contigs} \\
                   -p ./ \\
                   -o ${prefix}_comebin_bins \\
                   -n 8 \\
                   -t $task.cpus \\
                   -b 512 \\
                   -e 1024 \\
                   -c 1024
        fi

	"""
}
