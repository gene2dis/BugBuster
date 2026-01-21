process SOURMASH {
    container 'quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0'

    maxForks 24
    label 'process_single'

    publishDir "${params.output}/workflow/${meta.id}/sourmash_${db_name}", pattern: '*.with-lineages.csv'

    input:
        tuple val(meta), path(reads), path(db), path(tax_db)
	val db_name
        val sourmash_rank

    output:
        path("*.with-lineages.csv"), emit: sourmash_gather
	path("*_report.tsv"), emit: report

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        
        // Determine k-mer size based on taxonomic rank
        def kmer_map = [genus: 21, species: 31, strain: 51]
        def kmer = kmer_map[sourmash_rank.toLowerCase()] ?: 31
        
        // Determine input files based on single-end vs paired-end
        def input_reads = reads instanceof List && reads.size() > 2 && reads[2] != null ? 
            "${reads[0]} ${reads[1]} ${reads[2]}" : 
            reads instanceof List && reads.size() > 1 ? "${reads[0]} ${reads[1]}" : "${reads[0]}"

        """
        # Sketch reads with appropriate k-mer size
        sourmash sketch dna <(zcat ${input_reads}) -p k=${kmer},scaled=1000,abund -o ${prefix}.sig
        
        # Rename signature for clarity
        sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k${kmer}_renamed.sig
        
        # Gather matches against database (allow failure if no matches found)
        sourmash gather ${prefix}_k${kmer}_renamed.sig ${db} -k ${kmer} -o ${prefix}_smgather_${db_name}.csv || true
        
        # Check if gather produced results
        if [ -f "${prefix}_smgather_${db_name}.csv" ] && [ -s "${prefix}_smgather_${db_name}.csv" ]; then
            # Annotate with taxonomy using pre-prepared database
            sourmash tax annotate -g ${prefix}_smgather_${db_name}.csv -t ${tax_db}
            
            # Calculate classification statistics
            classified=\$(awk -F, '{sum += \$5} END {print sum}' ${prefix}_smgather_${db_name}.with-lineages.csv)
            unclassified=\$(echo "1 \${classified}" | awk '{print \$1 - \$2}')
        else
            # No matches found - create empty output files
            echo "query_filename,query_name,query_md5,query_bp,query_abundance,match_name,match_md5,gather_result_rank,f_match_orig,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,f_match,unique_intersect_bp,remaining_bp,query_containment_ani,match_containment_ani,average_containment_ani,max_containment_ani,n_unique_weighted_found,sum_weighted_found,total_weighted_hashes,lineage" > ${prefix}_smgather_${db_name}.with-lineages.csv
            classified=0
            unclassified=1
        fi
        
        # Generate report
        echo -e "Id\tSourmash DB\tUnclassified\tClassified" > ${meta.id}_${sourmash_rank}_${db_name}_report.tsv
        echo -e "${prefix}\t${db_name}\t\${unclassified}\t\${classified}" >> ${meta.id}_${sourmash_rank}_${db_name}_report.tsv
        
        # Cleanup temporary files
        rm -f *.sig
        """
}
