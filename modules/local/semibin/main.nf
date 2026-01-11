process SEMIBIN {
    container 'quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/semibin", pattern: '*semibin_bins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/raw_bins", mode: 'copy', pattern: '*semibin_output_bins'

    input:
        tuple val(meta), path(contigs), path(bam)
    
    output:
        tuple val(meta), path("*_semibin_output_bins"), emit: rosella

    script:
        def prefix = "${meta.id}"
        def env_model = "${params.semibin_env_model}"

        """
        min_length=1000
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' ${contigs})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs in ${contigs} are shorter than \${min_length}bp. Skipping SemiBin for sample ${prefix}." >&2
            mkdir -p ${prefix}_semibin_output_bins
            echo "Sample ${prefix}: All contigs < \${min_length}bp - SemiBin skipped" > ${prefix}_semibin_output_bins/SKIPPED.txt
        else
            SemiBin2 single_easy_bin \\
                     -i ${contigs} \\
                     -b ${bam} \\
                     -o ${prefix}_semibin_bins \\
                     --environment ${env_model} \\
                     --compression none \\
                     --threads $task.cpus

            mv ${prefix}_semibin_bins/output_bins ${prefix}_semibin_output_bins
            rm -rf ${prefix}_semibin_bins
        fi
        
        chmod 777 -R ${prefix}_semibin_output_bins
	"""
}

process SEMIBIN_COASSEMBLY {
    container 'quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0'

    label 'process_high'

    publishDir "${params.output}/workflow/co_assembly/semibin", pattern: 'co_assembly_semibin_bins'
    publishDir "${params.output}/contigs_and_bins/co_assembly/raw_bins", mode: 'copy', pattern: 'co_assembly_semibin_bins'

    input:
        path(bams_and_contigs)
    
    output:
        path("co_assembly_semibin_bins"), emit: semibin

    script:
        """
        contigs=`ls | grep -E '.+?filtered_contigs.+' | tr '\\n' ',' | sed 's/.\$//'`
        bam_list=`ls | grep -E '.+?all_reads.bam' | tr '\\n' ' ' | sed 's/.\$//'`

        min_length=1000
        contig_file=\$(echo \${contigs} | cut -d',' -f1)
        long_contigs=\$(awk -v min=\$min_length '/^>/ {if (seqlen >= min) count++; seqlen=0; next} {seqlen += length(\$0)} END {if (seqlen >= min) count++; print count+0}' \${contig_file})
        
        if [ "\$long_contigs" -eq 0 ]; then
            echo "WARNING: All contigs in coassembly are shorter than \${min_length}bp. Skipping SemiBin." >&2
            mkdir -p co_assembly_semibin_bins
            echo "Coassembly: All contigs < \${min_length}bp - SemiBin skipped" > co_assembly_semibin_bins/SKIPPED.txt
        else
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
            rm -rf output
        fi
        
        chmod 777 -R co_assembly_semibin_bins
	"""
}
