process PRODIGAL_BINS {
    container 'quay.io/biocontainers/prodigal:2.6.3--h031d066_8'

    label 'process_medium'

    publishDir "${params.output}/workflow/${meta.id}/prodigal_bins", pattern: '*_bins_proteins'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/refined_bins_proteins", mode: 'copy', pattern: '*_bins_proteins'

    input:
        tuple val(meta), path(metawrap)

    output:
        tuple val(meta), path("*_bins_proteins"), emit: prodigal_bins

    script:
        def prefix = "${meta.id}"

        """
        #!/bin/bash

        cp -rL ${metawrap} tmp_bins
	cd tmp_bins

        mkdir ${prefix}_bins_genes
        mkdir ${prefix}_bins_proteins

        N=$task.cpus
        i=1
        
	for file in *.fa; do 
             if (( \$i % \$N == 0 )); then
                wait
             fi
             ((i++))

             file_name=`echo \$file | sed 's/.fa//g'`
             prodigal -i \$file \\
                      -o ${prefix}_\${file_name}_genes.gff \\
                      -a ${prefix}_\${file_name}_proteins.faa \\
                      -p single &
        done; wait
        
        mv *_genes.gff ${prefix}_bins_genes
        mv *_proteins.faa ${prefix}_bins_proteins
	cd ..
        mv tmp_bins/${prefix}_bins_genes/ .
        mv tmp_bins/${prefix}_bins_proteins/ .

        rm -rf refined_bins/
	"""
}

process PRODIGAL_CONTIGS {
    container 'quay.io/biocontainers/prodigal:2.6.3--h031d066_8'

    label 'process_single'

    publishDir "${params.output}/workflow/${meta.id}/prodigal_contigs", pattern: '*contigs*'
    publishDir "${params.output}/contigs_and_bins/${meta.id}/contigs", mode: 'copy', pattern: '*contigs*'

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("*_contigs_proteins.faa"), emit: prodigal_contigs

    script:
        def prefix = "${meta.id}"

        """
        prodigal -i ${contigs} \\
                 -o ${prefix}_contigs_genes.gff \\
                 -a ${prefix}_contigs_proteins.faa \\
                 -p meta 
        """
}
