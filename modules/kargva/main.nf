process KARGVA {
    container 'quay.io/ffuentessantander/kargva:1.1'

    label 'process_low'

    input:
        tuple val(meta), path(reads)
        path(params.kargva_db)
    
    output:
        tuple val(meta), path("*_all_reads_KARGVA_mappedGenes.csv"), path(reads), emit: kargva_reads
        path("*_all_reads_KARGVA_mappedGenes.csv"), emit: kargva_reports

    script:
        def prefix = "${meta.id}"

        """
	if [[ ${reads[2]} != null ]]; then
		cat ${reads[0]} ${reads[1]} ${reads[2]} > ${prefix}_all_reads.fastq.gz
	else 
		cat ${reads[0]} ${reads[1]} > ${prefix}_all_reads.fastq.gz
	fi

	java -cp /bin/ KARGVA k:17 d:${params.kargva_db} -XX:ActiveProcessorCount=${task.cpus} -Xmx32GB ${prefix}_all_reads.fastq.gz

        if [[ ! -e ${prefix}_all_reads_KARGVA_mappedGenes.csv ]]; then
                echo "GeneIdx,KmerSNPHits,PercentGeneCovered,AverageKMerDepth" > ${prefix}_all_reads_KARGVA_mappedGenes.csv
                echo "NA,NA,NA,NA" >> ${prefix}_all_reads_KARGVA_mappedGenes.csv
        fi

        rm -f *KARGVA_mappedReads.csv
        rm -f ${prefix}_all_reads.fastq.gz
        """
}
