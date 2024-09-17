process FASTP {
    container 'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0'

    label 'process_low'

    input:
        tuple val(meta), path(reads)
    
    output:
        tuple val(meta), path("*_fastp.fastq.gz"), path("*_report.json"), emit: fastq

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_1_fastp.fastq.gz \\
            --out2 ${prefix}_2_fastp.fastq.gz \\
            -j ${prefix}_report.json \\
            $args \\
            --thread $task.cpus && \\

	if [[ ${reads[2]} != null ]]; then
		fastp \\
	            -i ${reads[2]} \\
	            -o ${prefix}_Singleton_fastp.fastq.gz \\
		    -j ${prefix}_Singleton_report.json \\
                    $args \\
		    --thread $task.cpus
	fi
        """
}
