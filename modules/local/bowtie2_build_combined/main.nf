/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOWTIE2 BUILD COMBINED INDEX Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Build a single Bowtie2 index from multiple FASTA files for efficient decontamination
----------------------------------------------------------------------------------------
*/

process BOWTIE2_BUILD_COMBINED {
    tag "${index_name}"
    container 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0'

    label 'process_medium'

    input:
    path(fasta_files)
    val index_name

    output:
    path("${index_name}_index"), emit: index
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Create output directory
    mkdir ${index_name}_index

    # Concatenate all FASTA files into a single file
    # Handle both compressed and uncompressed files
    for fasta in ${fasta_files}; do
        if [[ \$fasta == *.gz ]]; then
            zcat \$fasta >> combined_contaminants.fasta
        else
            cat \$fasta >> combined_contaminants.fasta
        fi
    done

    # Build Bowtie2 index from combined FASTA
    bowtie2-build \\
        --threads ${task.cpus} \\
        combined_contaminants.fasta \\
        ${index_name}_index/${index_name}

    # Clean up temporary combined FASTA
    rm combined_contaminants.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${index_name}_index
    touch ${index_name}_index/${index_name}.1.bt2
    touch ${index_name}_index/${index_name}.2.bt2
    touch ${index_name}_index/${index_name}.3.bt2
    touch ${index_name}_index/${index_name}.4.bt2
    touch ${index_name}_index/${index_name}.rev.1.bt2
    touch ${index_name}_index/${index_name}.rev.2.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.3
    END_VERSIONS
    """
}
