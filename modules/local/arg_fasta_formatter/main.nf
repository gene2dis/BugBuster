process ARG_FASTA_FORMATTER {
    tag "${meta.id}"
    label 'process_single'

    // The report processes' image (biopython is in it for this module); the
    // helper, extract_sequences.py, is in bin/
    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"

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
        python: 3.13.18
        biopython: 1.88
    END_VERSIONS
    """
}
