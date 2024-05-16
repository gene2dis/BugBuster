process AUTOMETA {
    container 'quay.io/biocontainers/autometa:2.2.0--pyh7cba7a3_0'

    cpus 30

    publishDir "${params.output}/workflow/co_assembly/autometa", pattern: 'autometa_bins'
    publishDir "${params.output}/Co_assembly/Bins", mode: 'copy', pattern: 'autometa_bins'

    input:
        path(bams_and_contigs)
	path(ncbi_db)
    
    output:
        path("autometa_bins"), emit: autometa

    script:
        def prefix = "co_assembly"
       
        """
        contigs=`ls | grep -E '.+?filtered_contigs.+' | tr '\\n' ',' | sed 's/.\$//'`
        bam_list=`ls | grep -E '.+?all_reads.bam' | tr '\\n' ' ' | sed 's/.\$//'`

        samtools merge ${prefix}_all_samples.bam \${bam_list}

	autometa-length-filter \\
	    --assembly \${contigs} \\
	    --cutoff 1000 \\
	    --output-fasta ${prefix}_contigs.filtered.fna \\
	    --output-stats ${prefix}_contigs.stats.tsv \\
	    --output-gc-content ${prefix}_contigs.gc_content.tsv

        autometa-coverage \\
            --assembly \${contigs} \\
            --bam ${prefix}_all_samples.bam \\
            --out ./${prefix}_contigs_coverages.tsv \\
            --cpus $task.cpus

	autometa-orfs \\
            --assembly ${prefix}_contigs.filtered.fna \\
            --output-nucls ./${prefix}_contigs_orfs.fna \\
            --output-prots ./${prefix}_contigs_orfs.faa \\
            --cpus $task.cpus

	autometa-markers \\
            --orfs ./${prefix}_contigs_orfs.faa \\
            --kingdom bacteria \\
            --hmmscan ./${prefix}_contigs.hmmscan.tsv \\
            --out ./${prefix}_contigs.markers.tsv \\
            --dbdir ${ncbi_db}/markers \\
            --hmmdb ${ncbi_db}/markers \\
            --cpus $task.cpus \\
            --seed 42

        diamond blastp \\
            --query ${prefix}_contigs_orfs.faa \\
            --db ${ncbi_db}/nr.dmnd \\
            --evalue 1e-5 \\
            --max-target-seqs 200 \\
            --threads $task.cpus \\
            --outfmt 6 \\
            --out ./${prefix}_contigs.blastp.tsv

        autometa-taxonomy-lca \\
            --blast ./${prefix}_contigs.blastp.tsv \\
            --dbdir ${ncbi_db} \\
            --lca-output ./${prefix}_contigs.lca.tsv \\
            --sseqid2taxid-output ./${prefix}_contigs.lca.sseqid2taxid.tsv \\
            --lca-error-taxids ./${prefix}_contigs.lca.errorTaxids.tsv

        autometa-taxonomy-majority-vote \\
            --lca ./${prefix}_contigs.lca.tsv \\
            --output ./${prefix}_contigs.votes.tsv \\
            --dbdir ${ncbi_db}

        autometa-taxonomy \\
            --votes ./${prefix}_contigs.votes.tsv \\
            --output ./ \\
            --assembly ./${prefix}_contigs.filtered.fna \\
            --prefix ${prefix}_contigs \\
            --split-rank-and-write superkingdom \\
	    --dbdir ${ncbi_db} \\
            --dbtype ncbi

        autometa-kmers \\
            --fasta ./${prefix}_contigs.bacteria.fna \\
            --kmers ./${prefix}_contigs.bacteria.kmers.tsv \\
            --size 5 \\
            --norm-method am_clr \\
            --norm-output ./${prefix}_contigs.bacteria.kmers.normalized.tsv \\
            --pca-dimensions 50 \\
            --embedding-method bhsne \\
            --embedding-output ./${prefix}_contigs.bacteria.kmers.embedded.tsv \\
            --cpus $task.cpus \\
            --seed 42

        autometa-binning \\
            --kmers ${prefix}_contigs.bacteria.kmers.embedded.tsv \\
            --coverages ${prefix}_contigs_coverages.tsv \\
            --gc-content ${prefix}_contigs.gc_content.tsv \\
            --markers ${prefix}_contigs.markers.tsv \\
            --output-binning ./${prefix}_contigs_binning.tsv \\
            --output-main ./${prefix}_contigs.main.tsv \\
            --clustering-method dbscan \\
            --completeness 20 \\
            --purity 90 \\
            --cov-stddev-limit 25 \\
            --gc-stddev-limit 5 \\
            --taxonomy ${prefix}_contigs.taxonomy.tsv \\
            --starting-rank superkingdom \\
            --rank-filter superkingdom \\
            --rank-name-filter bacteria \\
            --cpus $task.cpus

        autometa-binning-summary \\
            --binning-main ./${prefix}_contigs.main.tsv \\
            --markers ./${prefix}_contigs.markers.tsv \\
            --metagenome ./${prefix}_contigs.filtered.fna \\
            --output-stats ./${prefix}_contigs_summary.stats \\
            --output-taxonomy ./${prefix}_contigs_summary.tax \\
            --output-metabins autometa_bins

        cd autometa_bins/
	rm unclustered*
        for file in bin*.fna ; do 
            file_name=`echo \$file | sed 's/fna/fa/g'`
            cat \$file | sed -E '/>/ s/ .+//g' > \$file_name
            rm \$file
        done
	cd ../  
	"""
}

