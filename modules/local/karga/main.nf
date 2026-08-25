process KARGA {
    container 'quay.io/ffuentessantander/karga:1.1'

    label 'process_low'

    input:
        tuple val(meta), path(report), path(reads), path(karga_db)

    output:
        path("*_all_reads_KARGA_mappedGenes.csv"), emit: kargva
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        if [[ ${reads[2]} != null ]]; then
                cat ${reads[0]} ${reads[1]} ${reads[2]} > ${prefix}_all_reads.fastq.gz
        else
                cat ${reads[0]} ${reads[1]} > ${prefix}_all_reads.fastq.gz
        fi

        java -XX:ActiveProcessorCount=${task.cpus} -Xmx${task.memory.toGiga()}g -cp /bin/ KARGA k:17 d:${karga_db} r:n ${prefix}_all_reads.fastq.gz
        rm -f ${prefix}_all_reads.fastq.gz

        if [[ ! -e ${prefix}_all_reads_KARGA_mappedGenes.csv ]]; then
                echo "GeneIdx,PercentGeneCovered,AverageKMerDepth" > ${prefix}_all_reads_KARGA_mappedGenes.csv
                echo "NA,NA,NA" >> ${prefix}_all_reads_KARGA_mappedGenes.csv
        fi

        # KARGA is a bare Java class with no version flag: version from the
        # container tag
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            karga: 1.1
            java: \$(java -version 2>&1 | head -n 1 | sed 's/.*"\\(.*\\)".*/\\1/')
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"

        """
        echo "GeneIdx,PercentGeneCovered,AverageKMerDepth" > ${prefix}_all_reads_KARGA_mappedGenes.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            karga: 1.1
            java: 18.0.2.1
        END_VERSIONS
        """
}
