process COMPUTE_KMER_FREQ {
    container 'quay.io/ffuentessantander/kmer_freq:1.0'

    cpus 1

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("*_kmer_freq.tsv"), emit: kmer_freq

    script:
        def prefix = "${meta.id}"

        """
        calc.kmerfreq.pl \\
                  -i ${contigs} \\
                  -o ${prefix}_kmer_freq.tsv
	"""
}
