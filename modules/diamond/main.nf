process PRODIGAL {
    container 'quay.io/biocontainers/prodigal:2.6.3--h031d066_8'

    cpus 10

    publishDir "${params.output}/workflow/${meta.id}/prodigal", mode: 'copy', pattern: '*_bins_proteins'

    input:
        tuple val(meta), path(metawrap)

    output:
        path("*_bins_proteins"), emit: gtdb_tk

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
