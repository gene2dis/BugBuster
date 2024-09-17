process ARG_FASTA_FORMATTER {

    container 'quay.io/ffuentessantander/biopython_utils:1.0'

    // publishDir "${params.output}/Arg_prediction/Contig_level", mode: 'copy', pattern: 'Contig_tax_and_arg_prediction.tsv'

    cpus 1

    input:
        path(arg_bin_data)

    output:
        path("*_deeparg_mapped.faa"), emit: arg_reports

    script: 

    """
    cp -rL $carpeta1/* . 
    cp -rL $carpeta2/* .
    for file in *_deep_arg.out.mapping.ARG; do \\
             file_name=`echo \$file | sed 's/_deep_arg.out.mapping.ARG//g'`
             protein_file=`echo \$file | sed 's/_deep_arg.out.mapping.ARG/_proteins.faa/g'`
             cat \$file | tail -n +2 | awk '{print \$4}' > tmp_list  
             if [[ `cat tmp_list` != "" ]]; then python extract_sequences.py \$protein_file tmp_list > \${file_name}_deeparg_mapped.faa; fi
             rm -f tmp_list
    done
    rm -f *mapping.potential.ARG
    rm -f *out.mapping.ARG
    rm -f *_proteins.faa
    """

}
