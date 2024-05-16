process SEMIBIN {
    container 'quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0'

    cpus 10

    publishDir "${params.output}/workflow/co_assembly/semibin", pattern: 'co_assembly_semibin_bins'
    publishDir "${params.output}/Co_assembly/Bins", mode: 'copy', pattern: 'co_assembly_semibin_bins'

    input:
        path(bams_and_contigs)
    
    output:
        path("co_assembly_semibin_bins"), emit: semibin

    script:
        """
        contigs=`ls | grep -E '.+?filtered_contigs.+' | tr '\\n' ',' | sed 's/.\$//'`
        bam_list=`ls | grep -E '.+?all_reads.bam' | tr '\\n' ' ' | sed 's/.\$//'`

        SemiBin2 generate_sequence_features_single \\
                 -i \${contigs} \\
                 -b \${bam_list} \\
                 -o contig_output \\
                 --threads $task.cpus

        SemiBin2 train_self \\
                 --data contig_output/data.csv \\
                 --data-split contig_output/data_split.csv \\
                 -o contig_output \\
                 --threads $task.cpus

        SemiBin2 bin_short \\
                 -i \${contigs} \\
                 --model contig_output/model.h5 \\
                 --data contig_output/data.csv \\
                 -o output \\
                 --compression none \\
                 --threads $task.cpus

        mv output/output_bins co_assembly_semibin_bins
        chmod 777 -R co_assembly_semibin_bins
	"""
}
