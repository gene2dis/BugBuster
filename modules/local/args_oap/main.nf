process ARGS_OAP {
    tag "${meta.id}"
    // process_medium (8 CPUs, 24 h): stage one on 55-70 M read pairs took 6-8 h with process_low's 4 CPUs,
    // and one sample hit its 8 h limit
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/args_oap:4.0.0--pyhdfd78af_0' :
        'quay.io/biocontainers/args_oap:4.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_args_oap"), emit: args_oap
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix  = "${meta.id}"
    def singles = reads instanceof List && reads.size() > 2 ? reads[2] : null
    """
    # args_oap globs its input folder, so the reads get a folder of their own;
    # links instead of copies (it reads through them, performance review E6).
    # The singleton file is a sample of its own for ARGs-OAP (<id>_tmp_reads_S);
    # ARGS_OAP_REPORT adds it to the pairs
    mkdir tmp_reads
    ln -s ../${reads[0]} tmp_reads/${prefix}_tmp_reads_R1.fastq.gz
    ln -s ../${reads[1]} tmp_reads/${prefix}_tmp_reads_R2.fastq.gz
    ${singles ? "ln -s ../${singles} tmp_reads/${prefix}_tmp_reads_S.fastq.gz" : ''}

    # The image ships SARG, KO30 and GG85 as FASTA and stage one indexes them
    # inside its own package directory on first use, which a read-only image
    # (singularity/apptainer) refuses (design doc Q21). A copy of the package
    # (~35 MB) in the task directory is run instead and indexed there (~12 s);
    # stock code, so nothing to guard against an upstream update
    pkg=\$(python -c 'import args_oap, os; print(os.path.dirname(args_oap.__file__))')
    mkdir pylib
    cp -rL "\$pkg" pylib/
    export PYTHONPATH="\$PWD/pylib\${PYTHONPATH:+:\$PYTHONPATH}"
    case \$(python -c 'import args_oap; print(args_oap.__file__)') in
        "\$PWD"/pylib/*) ;;
        *) echo "ERROR: args_oap is not imported from the task directory's copy" >&2; exit 1 ;;
    esac

    args_oap stage_one \\
        -i tmp_reads \\
        -o ${prefix}_args_oap \\
        -f fastq.gz \\
        -t ${task.cpus} \\
        ${args} > stage_one.log 2>&1 || { cat stage_one.log >&2; exit 1; }
    cat stage_one.log >&2
    # stage one reports a failure on an input file as a warning and drops the file
    if grep -q 'Something is wrong with' stage_one.log; then
        echo "ERROR: ARGs-OAP stage one dropped an input file of ${prefix} (see the log above)" >&2
        exit 1
    fi

    # Stage two normalises by the metadata rows with reads, 16S and cells all
    # above 0; without one the sample keeps metadata.txt only. Its outcome goes
    # to stage_two_status.txt for ARGS_OAP_REPORT: ok, no_sarg_hit (no read
    # left after the SARG search: stage two exits 2) or not_normalised
    if awk -F'\t' 'NR > 1 && \$2 > 0 && \$3 > 0 && \$4 > 0 { found = 1 } END { exit !found }' ${prefix}_args_oap/metadata.txt; then
        status=0
        args_oap stage_two \\
            -i ${prefix}_args_oap \\
            -t ${task.cpus} \\
            ${args2} > stage_two.log 2>&1 || status=\$?
        cat stage_two.log >&2
        if [ "\$status" -eq 0 ]; then
            echo ok > ${prefix}_args_oap/stage_two_status.txt
        elif [ "\$status" -eq 2 ] && grep -qE 'No target sequence (detected|remained)' stage_two.log; then
            echo no_sarg_hit > ${prefix}_args_oap/stage_two_status.txt
        else
            exit "\$status"
        fi
    else
        echo "ARGs-OAP: no 16S or cell estimate for ${prefix}, stage two skipped" >&2
        echo not_normalised > ${prefix}_args_oap/stage_two_status.txt
    fi

    rm -rf tmp_reads pylib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        args_oap: \$(args_oap --version 2>&1)
        diamond: \$(diamond --version 2>&1 | sed 's/diamond version //')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    mkdir ${prefix}_args_oap
    printf 'sample\\tnRead\\tn16S\\tnCell\\n${prefix}_tmp_reads\\t2000\\t2.0\\t1.0\\n' > ${prefix}_args_oap/metadata.txt
    printf 'qseqid\\tsseqid\\tpident\\tlength\\tqlen\\tslen\\tevalue\\tbitscore\\tqcov\\tsample\\tscov\\tgene\\tsubtype\\ttype\\tcount\\tcopy\\trpk\\trpksum\\n' > ${prefix}_args_oap/blastout.filtered.txt
    touch ${prefix}_args_oap/extracted.fa ${prefix}_args_oap/normalized_cell.type.txt
    echo no_sarg_hit > ${prefix}_args_oap/stage_two_status.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        args_oap: 4.0.0
        diamond: 2.1.13
    END_VERSIONS
    """
}
