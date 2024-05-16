process BEST_HITS_TO_FASTA {
    container 'quay.io/ffuentessantander/r_reports:1.1'

    input:
        path(mmseq_data)

    output:
        path("All_mmseq_results.faa"), emit: faa

    script:

        """
        Rscript /mnt/Best_hits_search.R 
        cat *.faa > All_mmseq_results.faa 
	"""
}
