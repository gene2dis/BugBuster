process KARGVA {
    container 'quay.io/ffuentessantander/kargva:1.0'

    cpus 1

    publishDir "${params.output}/workflow/${meta.id}/kargva", pattern: '*KARGVA_mappedGenes.csv'

    input:
        tuple val(meta), path(reads)
    
    output:
        tuple val(meta), path("*_all_reads_KARGVA_mappedGenes.csv"), path("*all_reads.fastq.gz"), emit: kargva_reads
        path("*_all_reads_KARGVA_mappedGenes.csv"), emit: kargva_reports

    script:
        def prefix = "${meta.id}"

        """
	mv /mnt/* .
	if [[ ${reads[2]} != null ]]; then
		cat ${reads[0]} ${reads[1]} ${reads[2]} > ${prefix}_all_reads.fastq.gz
	else 
		cat ${reads[0]} ${reads[1]} > ${prefix}_all_reads.fastq.gz
	fi
	java KARGVA k:17 -Xmx32GB ${prefix}_all_reads.fastq.gz
        rm -f *KARGVA_mappedReads.csv
        """
}
