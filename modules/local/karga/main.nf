process KARGA {
    container 'quay.io/ffuentessantander/karga:1.1'

    label 'process_low'

    input:
        tuple val(meta), path(report), path(reads), path(karga_db)

    output:
        path("*_all_reads_KARGA_mappedGenes.csv"), emit: kargva

    script:
        def prefix = "${meta.id}"

        """
        if [[ ${reads[2]} != null ]]; then
                cat ${reads[0]} ${reads[1]} ${reads[2]} > ${prefix}_all_reads.fastq.gz
        else
                cat ${reads[0]} ${reads[1]} > ${prefix}_all_reads.fastq.gz
        fi

        java -cp /bin/ KARGA k:17 d:${karga_db} -XX:ActiveProcessorCount=${task.cpus} r:n -Xmx32GB ${prefix}_all_reads.fastq.gz
        rm -f ${prefix}_all_reads.fastq.gz

        if [[ ! -e ${prefix}_all_reads_KARGA_mappedGenes.csv ]]; then
                echo "GeneIdx,PercentGeneCovered,AverageKMerDepth" > ${prefix}_all_reads_KARGA_mappedGenes.csv
                echo "NA,NA,NA" >> ${prefix}_all_reads_KARGA_mappedGenes.csv
        fi
        """
}
