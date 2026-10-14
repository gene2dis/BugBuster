/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWAMEM3 BUILD COMBINED INDEX Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Build a single bwa-mem3 index from multiple FASTA files (phiX + host) for
    single-pass decontamination (--read_aligner bwamem3)

    Peak memory is ~12 bytes per base of doubled text (~75 GB for CHM13 + phiX).
    config/modules.config sizes the request from the staged FASTA sizes,
    counting a gzipped file at 4x (the input is usually gzipped, and its size
    alone would undershoot by the compression ratio). --max-memory is pinned to the task
    allocation because bwa-mem3 otherwise budgets against the whole node's RAM
    (it is only cgroup-aware); a budget that cannot fit the reference makes
    bwa-mem3 exit 3 before building, which config/modules.config retries.
----------------------------------------------------------------------------------------
*/

process BWAMEM3_BUILD_COMBINED {
    tag "${index_name}"
    container 'community.wave.seqera.io/library/bwa-mem3_samtools:65dfdfc4140db421'

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
    def args = task.ext.args ?: ''
    def max_memory = task.memory ? "--max-memory ${task.memory.toMega()}M" : ''
    """
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

    # Exit 3: the reference does not fit in the memory budget (checked before
    # building). The retry raises it until --max_memory caps it
    status=0
    bwa-mem3 index \\
        ${args} \\
        -t ${task.cpus} \\
        ${max_memory} \\
        -p ${index_name}_index/${index_name} \\
        combined_contaminants.fasta || status=\$?
    if [ \$status -eq 3 ]; then
        echo "ERROR: building the bwa-mem3 index needs more memory than ${task.memory} (a human host needs ~75 GB). Raise --max_memory, pass a pre-built bwa-mem3 index with --custom_decontamination_index, or use --read_aligner bowtie2." >&2
    fi
    [ \$status -eq 0 ] || exit \$status

    rm combined_contaminants.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: \$(bwa-mem3 version | sed -nE '1 s/^([0-9]+(\\.[0-9]+)+).*/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${index_name}_index
    touch ${index_name}_index/${index_name}.amb
    touch ${index_name}_index/${index_name}.ann
    touch ${index_name}_index/${index_name}.bwt.2bit.64
    touch ${index_name}_index/${index_name}.pac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem3: 0.14.0
    END_VERSIONS
    """
}
