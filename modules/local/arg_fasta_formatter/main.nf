process ARG_FASTA_FORMATTER {

    container 'quay.io/ffuentessantander/biopython_utils:1.0'

    label 'process_single'

    input:
        tuple val(meta), path(proteins), path(deeparg)

    output:
        path("*_ARGs.faa"), emit: arg_reports
        path "versions.yml", emit: versions

    script:

    """
    cp -r ${proteins}/* . 
    cp -r ${deeparg}/* .
    for file in *_deep_arg.out.mapping.ARG; do \\
             file_name=`echo \$file | sed 's/_deep_arg.out.mapping.ARG//g'`
             protein_file=`echo \$file | sed 's/_deep_arg.out.mapping.ARG/_proteins.faa/g'`
             cat \$file | tail -n +2 | awk '{print \$4"\t"\$1"\t"\$5"\t"\$6}' > tmp_list  
             if [[ `cat tmp_list` != "" ]]; then extract_sequences.py \$protein_file tmp_list; fi
             rm -f tmp_list
    done
    rm -f *mapping.potential.ARG
    rm -f *out.mapping.ARG
    rm -f *_proteins.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    // Per-sample name: a fixed name collides when CLUSTERING collects the
    // outputs of every sample into one task
    """
    touch ${meta.id}_stub_ARGs.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.8.10
        biopython: 1.78
    END_VERSIONS
    """
}
