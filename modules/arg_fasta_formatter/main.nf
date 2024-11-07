process ARG_FASTA_FORMATTER {

    container 'quay.io/ffuentessantander/biopython_utils:1.0'

    label 'process_single'

    input:
        tuple val(meta), path(proteins), path(deeparg)

    output:
        path("*_ARGs.faa"), emit: arg_reports

    script: 

    """
    cp -rL ${proteins}/* . 
    cp -rL ${deeparg}/* .
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
    """
}
