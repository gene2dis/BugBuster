process BBMAP {
    container 'quay.io/biocontainers/bbmap:39.06--h92535d8_0'

    label 'process_single'

    publishDir "${params.output}/workflow/${meta.id}/bbmap", pattern: '*_filtered_contigs.fa'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/contigs", mode: 'copy', pattern: '*_filtered_contigs.fa'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/contigs", mode: 'copy', pattern: '*_contig.stats'

    input:
        tuple val(meta), path(reads), path(contigs)

    output:
        tuple val(meta), path(reads), path("*_filtered_contigs.fa"), emit: bbmap
        tuple val(meta), path("*_filtered_contigs.fa"), emit: bbmap_contigs
        path("*_contig.stats"), emit: reports

    script:
        def prefix = "${meta.id}"
        def contig_length = "${params.bbmap_lenght}"

        """
	reformat.sh in=${contigs} out=${prefix}_filtered_contigs.fa minlength=${contig_length}	
	stats.sh in=${contigs} out=${prefix}_contig.stats
	"""
}

process BBMAP_COASSEMBLY {
    container 'quay.io/biocontainers/bbmap:39.06--h92535d8_0'

    label 'process_single'

    publishDir "${params.output}/workflow/co_assembly/bbmap", pattern: '*_filtered_contigs.fa'
    publishDir "${params.output}/contigs_and_bins/co_assembly/contigs", mode: 'copy', pattern: '*_contig.stats'
    publishDir "${params.output}/contigs_and_bins/co_assembly/contigs", mode: 'copy', pattern: '*_filtered_contigs.fa'

    input:
        tuple path(reads), path(contigs)

    output:
        path("*_filtered_contigs.fa"), emit: contigs
        path("*_contig.stats"), emit: reports

    script:
        def contig_length = "${params.bbmap_lenght}"

        """
	reformat.sh in=${contigs} out=co_assembly_filtered_contigs.fa minlength=${contig_length}	
	stats.sh in=${contigs} out=co_assembly_contig.stats
	"""
}
