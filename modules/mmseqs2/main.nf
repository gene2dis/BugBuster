process MMSEQS2 {
    container 'quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_1'
        
    publishDir = [path: { "${params.output}/workflow/${meta.id}/mmseqs2" }, pattern: "*_mmseqs_results"]

    cpus 10

    input:
        tuple val(meta), path(prodigal)
        path(arg_db)

    output:
        path("*_mmseqs_results"), emit: mmseqs2

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        cp -rL ${prodigal} tmp_bins
        mkdir ${prefix}_mmseqs_results
        cd tmp_bins
        for bin in *.faa; do
               bin_name=`echo \$bin | sed -E 's/_proteins.faa//g'` 

               mmseqs easy-search ../${arg_db} \\
                           \${bin} \\
                           \${bin_name}_mmseq.m8 \\
                           tmp \\
                           --threads $task.cpus \\
                           $args

              cat \${bin_name}_mmseq.m8 | sed -E "s/-\t/\$bin_name\t/g" | sed -E "s/empty/sample_bin/g" | grep -v 'RequiresSNPConfirmation' > \${bin_name}_mmseqs_result.tsv
              rm \${bin_name}_mmseq.m8
        done
        mv *_mmseqs_result.tsv ../${prefix}_mmseqs_results/
        cd ..
        rm -rf tmp_bins
	"""
}
