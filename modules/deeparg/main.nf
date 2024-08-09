process DEEPARG_BINS {
    container 'quay.io/biocontainers/deeparg:1.0.4--pyhdfd78af_0'
        
    publishDir = [path: { "${params.output}/workflow/${meta.id}/deeparg_bins" }, pattern: "*_deeparg_results"]
    publishDir = [path: { "${params.output}/ARG_prediction/Bin_level/${meta.id}/DeepARG" }, mode: 'copy', pattern: "*_deeparg_results"]

    label 'process_low'

    input:
        tuple val(meta), path(prodigal_bins)
        path(deeparg_db)

    output:
        tuple val(meta), path("*_deeparg_results"), emit: deeparg_bins

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        cat /usr/local/lib/python2.7/site-packages/deeparg/entry.py | \\
        sed -e "s/                        ])/                        \\'-p $task.cpus\\'])/g" \\
        -e "s/output_file...align'/output_file\\+\\'.align\\',/g" > entry.py
        rm -f /usr/local/lib/python2.7/site-packages/deeparg/entry.py
        mv entry.py /usr/local/lib/python2.7/site-packages/deeparg/

        cp -rL ${prodigal_bins} tmp_bins
        mkdir ${prefix}_deeparg_results
        cd tmp_bins
        for bin_prot in *.faa; do
               bin_name=`echo \$bin_prot | sed -E 's/_proteins.faa//g'` 

               deeparg predict \\
                          -d ../${deeparg_db} \\
                          --model LS \\
                          --type prot \\
                          $args \\
                          --input \${bin_prot} \\
                          --out \${bin_name}_deep_arg.out

        done
        mv *.mapping.ARG ../${prefix}_deeparg_results/
        cd ..
        rm -rf tmp_bins
	"""
}

process DEEPARG_CONTIGS {
    container 'quay.io/biocontainers/deeparg:1.0.4--pyhdfd78af_0'

    publishDir = [path: { "${params.output}/workflow/${meta.id}/deeparg_contigs" }, pattern: "*_deeparg_results"]
    publishDir = [path: { "${params.output}/ARG_prediction/Contig_level/${meta.id}/DeepARG" }, mode: 'copy', pattern: "*_deeparg_results"]

    label 'process_low'

    input:
        tuple val(meta), path(prodigal_contigs)
        path(deeparg_db)

    output:
        tuple val(meta), path("*.ARG"), emit: deeparg
        path("*.ARG"), emit: only_deeparg

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        cat /usr/local/lib/python2.7/site-packages/deeparg/entry.py | sed -e "s/                        ])/                        \\'-p $task.cpus\\'])/g" -e "s/output_file...align'/output_file\\+\\'.align\\',/g" > entry.py
        rm -f /usr/local/lib/python2.7/site-packages/deeparg/entry.py
        mv entry.py /usr/local/lib/python2.7/site-packages/deeparg/

        deeparg predict \\
                -d ${deeparg_db} \\
                --model LS \\
                --type prot \\
                $args \\
                --input ${prodigal_contigs} \\
                --out ${prefix}_contigs_deep_arg.out
        """
}
