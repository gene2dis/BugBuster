process ARGS_OAP {

    container 'quay.io/ffuentessantander/args_oap:3.2.4'

    label 'process_low'

    input:
        tuple val(meta), path(reads)

    output:
        path("*args_oap_s1_out"), emit: args_oap_s1
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

    """
    mkdir tmp_reads
    cp -rL ${reads[0]} tmp_reads/${prefix}_tmp_reads_R1.fastq.gz
    cp -rL ${reads[1]} tmp_reads/${prefix}_tmp_reads_R2.fastq.gz

    if [ ${reads[2]} != null ]; then
        cp -rL ${reads[2]} tmp_reads/${prefix}_tmp_reads_S.fastq.gz
    fi

    args_oap stage_one -i tmp_reads \\
                       -o ${prefix}_args_oap_s1_out \\
                       -f fastq.gz \\
                       -t $task.cpus

    rm -rf tmp_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        args_oap: \$(pip show args_oap 2>/dev/null | sed -n 's/Version: //p' || echo 3.2.4)
    END_VERSIONS
    """

    stub:
        def prefix = "${meta.id}"

    """
    mkdir ${prefix}_args_oap_s1_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        args_oap: 3.2.4
    END_VERSIONS
    """
}
