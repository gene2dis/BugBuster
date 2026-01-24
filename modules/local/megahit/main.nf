/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MEGAHIT Assembly Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenome assembly using MEGAHIT
    
    Supports both per-sample assembly and co-assembly modes:
    - Per-sample: Assembles reads from individual samples
    - Co-assembly: Assembles pooled reads from multiple samples
    
    Input:
        tuple val(meta), path(reads)
    
    Output:
        contigs_and_reads: tuple val(meta), path(reads), path(contigs)
        contigs_only: tuple val(meta), path(contigs)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/megahit:1.2.9--h43eeafb_4' :
        'quay.io/biocontainers/megahit:1.2.9--h43eeafb_4' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), path("${meta.id}_contigs.fa"), emit: contigs_and_reads
    tuple val(meta), path("${meta.id}_contigs.fa")            , emit: contigs_only
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    
    // Input validation
    if (!meta || !meta.id) {
        error "MEGAHIT: meta.id is required"
    }
    if (!reads || reads.size() == 0) {
        error "MEGAHIT: No read files provided for sample ${meta.id}"
    }
    
    // Determine if this is coassembly mode based on meta.id
    def is_coassembly = (meta.id == 'coassembly')
    
    if (is_coassembly) {
        // Co-assembly mode: reads are already collected, need to separate R1/R2/singletons
        """
        set -euo pipefail
        
        # Separate reads by type for co-assembly
        R1_files=()
        R2_files=()
        Singleton_files=()
        
        for read_file in ${reads}; do
            if [[ \$read_file =~ _R1_|_1\\. ]]; then
                R1_files+=("\$read_file")
            elif [[ \$read_file =~ _R2_|_2\\. ]]; then
                R2_files+=("\$read_file")
            elif [[ \$read_file =~ [Ss]ingleton|[Uu]npaired ]]; then
                Singleton_files+=("\$read_file")
            fi
        done
        
        # Build comma-separated lists
        R1_list=\$(IFS=,; echo "\${R1_files[*]}")
        R2_list=\$(IFS=,; echo "\${R2_files[*]}")
        Singleton_list=\$(IFS=,; echo "\${Singleton_files[*]}")
        
        # Validate that we have R1 and R2 files
        if [ \${#R1_files[@]} -eq 0 ] || [ \${#R2_files[@]} -eq 0 ]; then
            echo "ERROR: Could not identify R1 and R2 files for co-assembly" >&2
            echo "Files provided: ${reads}" >&2
            exit 1
        fi
        
        # Run MEGAHIT with appropriate arguments
        if [ -n "\$Singleton_list" ]; then
            megahit \\
                -1 \$R1_list \\
                -2 \$R2_list \\
                --read \$Singleton_list \\
                -t ${task.cpus} \\
                ${args} \\
                -o ${prefix}_assembly
        else
            megahit \\
                -1 \$R1_list \\
                -2 \$R2_list \\
                -t ${task.cpus} \\
                ${args} \\
                -o ${prefix}_assembly
        fi
        
        mv ${prefix}_assembly/final.contigs.fa ${prefix}_contigs.fa
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(megahit --version 2>&1 | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    } else {
        // Per-sample assembly mode
        """
        set -euo pipefail
        
        R1=${reads[0]}
        R2=${reads[1]}
        
        # Validate input files exist
        if [ ! -f "\$R1" ] || [ ! -f "\$R2" ]; then
            echo "ERROR: Required read files not found" >&2
            echo "R1: \$R1" >&2
            echo "R2: \$R2" >&2
            exit 1
        fi
        
        # Check for singleton reads
        if [ ${reads.size()} -gt 2 ]; then
            SINGLETON=${reads[2]}
            megahit \\
                -1 \$R1 \\
                -2 \$R2 \\
                --read \$SINGLETON \\
                -t ${task.cpus} \\
                ${args} \\
                -o ${prefix}_assembly
        else
            megahit \\
                -1 \$R1 \\
                -2 \$R2 \\
                -t ${task.cpus} \\
                ${args} \\
                -o ${prefix}_assembly
        fi
        
        mv ${prefix}_assembly/final.contigs.fa ${prefix}_contigs.fa
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(megahit --version 2>&1 | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }

    stub:
    prefix = meta.id
    """
    # Create a minimal test contig file
    cat > ${prefix}_contigs.fa << 'EOF'
>contig_1 length=1234
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>contig_2 length=5678
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: 1.2.9
    END_VERSIONS
    """
}
