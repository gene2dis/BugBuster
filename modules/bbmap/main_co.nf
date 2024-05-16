process BBMAP {
    container 'quay.io/biocontainers/bbmap:39.06--h92535d8_0'

    cpus 5

    publishDir "${params.output}/workflow/co_assembly/bbmap", pattern: '*_filtered_contigs.fa'
    publishDir "${params.output}/Co_assembly/contigs", mode: 'copy', pattern: '*_contig.stats'
    publishDir "${params.output}/Co_assembly/contigs", mode: 'copy', pattern: '*_filtered_contigs.fa'

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
