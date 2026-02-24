process KARGA {
    container 'quay.io/ffuentessantander/karga:1.1'
    containerOptions '-v /bin/ps:/usr/bin/ps:ro -v /bin/ps:/bin/ps:ro -v /lib/x86_64-linux-gnu/libprocps.so.8:/lib64/libprocps.so.8:ro'

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

        java -XX:ActiveProcessorCount=${task.cpus} -Xmx${task.memory.toGiga()}g -cp /bin/ KARGA k:17 d:${karga_db} r:n ${prefix}_all_reads.fastq.gz
        rm -f ${prefix}_all_reads.fastq.gz

        if [[ ! -e ${prefix}_all_reads_KARGA_mappedGenes.csv ]]; then
                echo "GeneIdx,PercentGeneCovered,AverageKMerDepth" > ${prefix}_all_reads_KARGA_mappedGenes.csv
                echo "NA,NA,NA" >> ${prefix}_all_reads_KARGA_mappedGenes.csv
        fi
        """
}
