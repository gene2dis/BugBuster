/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOWTIE2_SAMTOOLS Alignment Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads to assembled contigs using Bowtie2 and generate sorted BAM files
    
    One sample per task, against its own assembly or the co-assembly
    Handles paired-end reads with optional singleton/unpaired reads
    Uses streaming SAM→BAM conversion to minimize I/O
    
    Input:
        tuple val(meta), path(reads), path(contigs)
    
    Output:
        contigs_and_bam: tuple val(meta), path(contigs), path(bam)
        bam_only: tuple val(meta), path(bam)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process BOWTIE2_SAMTOOLS {
    tag "${meta.id}"
    label 'process_high'
    
    // Seqera Containers bowtie2 2.5.5 + samtools 1.24, pinned for every engine:
    // its singularity build has no gzip, zcat or find (CONTRIBUTING)
    container 'community.wave.seqera.io/library/bowtie2_samtools:6243cc3671a7ee48'

    input:
    tuple val(meta), path(reads), path(contigs)

    output:
    tuple val(meta), path(contigs), path("${meta.id}_all_reads.bam"), emit: contigs_and_bam
    tuple val(meta), path("${meta.id}_all_reads.bam")              , emit: bam_only
    path "${meta.id}_empty_contig_report.txt", optional: true      , emit: empty_contig_report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    
    // Input validation
    if (!meta || !meta.id) {
        error "BOWTIE2_SAMTOOLS: meta.id is required"
    }
    if (!reads || reads.size() == 0) {
        error "BOWTIE2_SAMTOOLS: No read files provided for sample ${meta.id}"
    }
    if (!contigs) {
        error "BOWTIE2_SAMTOOLS: No contig file provided for sample ${meta.id}"
    }
    
    // One sample per task (R1, R2[, singleton] by position). Under
    // co-assembly each sample is mapped on its own against the co-assembly
    // and the pooled BAM is a samtools merge of these (issue #28).
    """
    set -euo pipefail
    
    # Validate input files
    if [ ! -f "${reads[0]}" ] || [ ! -f "${reads[1]}" ]; then
        echo "ERROR: Required read files not found" >&2
        echo "R1: ${reads[0]}" >&2
        echo "R2: ${reads[1]}" >&2
        exit 1
    fi
    
    if [ ! -f "${contigs}" ]; then
        echo "ERROR: Contig file not found: ${contigs}" >&2
        exit 1
    fi
    
    # Check if contig file is empty (no sequences after filtering)
    if ! grep -q '^>' "${contigs}"; then
        echo "WARNING: Contig file is empty for sample ${prefix} - creating empty BAM" >&2
        
        # Create minimal empty BAM with header only
        printf '@HD\\tVN:1.6\\tSO:coordinate\\n' > empty.sam
        samtools view -bS empty.sam > ${prefix}_all_reads.bam
        
        # Generate empty contig tracking report
        printf 'sample_id\\t${prefix}\\n' > ${prefix}_empty_contig_report.txt
        printf 'stage\\tBOWTIE2_SAMTOOLS\\n' >> ${prefix}_empty_contig_report.txt
        printf 'empty_contig_detected\\tTRUE\\n' >> ${prefix}_empty_contig_report.txt
        printf 'timestamp\\t%s\\n' "\$(date -u +"%Y-%m-%d %H:%M:%S UTC")" >> ${prefix}_empty_contig_report.txt
        printf 'action_taken\\tCreated empty BAM file\\n' >> ${prefix}_empty_contig_report.txt
        
        # Generate versions file
        printf '"%s":\\n' "${task.process}" > versions.yml
        printf '    bowtie2: %s\\n' "\$(bowtie2 --version 2>/dev/null | head -n 1 | sed 's/.*version //; s/ .*//')" >> versions.yml
        printf '    samtools: %s\\n' "\$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')" >> versions.yml
        
        exit 0
    fi
    
    # Build Bowtie2 index
    bowtie2-build \\
        --threads ${task.cpus} \\
        ${contigs} \\
        ${prefix}_index
    
    # Align paired-end reads and convert to sorted BAM via pipe
    bowtie2 \\
        ${args} \\
        -x ${prefix}_index \\
        -p ${task.cpus} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        2> ${prefix}_paired.log \\
        | samtools sort \\
            ${args2} \\
            -@ ${task.cpus} \\
            -o ${prefix}_paired.bam -
    
    # Handle singleton reads if present
    if [ ${reads.size()} -gt 2 ]; then
        bowtie2 \\
            ${args} \\
            -x ${prefix}_index \\
            -p ${task.cpus} \\
            -U ${reads[2]} \\
            2> ${prefix}_singleton.log \\
            | samtools sort \\
                ${args2} \\
                -@ ${task.cpus} \\
                -o ${prefix}_singleton.bam -
        
        # CLEANUP: Remove Bowtie2 index files after singleton alignment
        rm -f ${prefix}_index*.bt2 2>/dev/null || true
        
        # Merge paired and singleton BAMs
        samtools merge \\
            -@ ${task.cpus} \\
            ${prefix}_all_reads.bam \\
            ${prefix}_paired.bam \\
            ${prefix}_singleton.bam
        
        # CLEANUP: Remove intermediate BAM files after merging
        rm -f ${prefix}_paired.bam ${prefix}_singleton.bam 2>/dev/null || true
    else
        # No singletons - rename paired BAM
        mv ${prefix}_paired.bam ${prefix}_all_reads.bam
        # CLEANUP: Remove Bowtie2 index files
        rm -f ${prefix}_index*.bt2 2>/dev/null || true
    fi
    
    # CLEANUP: Remove log files
    rm -f ${prefix}_paired.log ${prefix}_singleton.log 2>/dev/null || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>/dev/null | head -n 1 | sed 's/.*version //; s/ .*//')
        samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    prefix = meta.id
    """
    # Create a minimal valid BAM file header
    # Note: In real stub testing, you might want to use samtools to create a proper empty BAM
    echo -e "@HD\tVN:1.6\tSO:coordinate" > header.sam
    echo -e "@SQ\tSN:contig_1\tLN:1000" >> header.sam
    echo -e "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.5.5" >> header.sam
    
    # Convert to BAM (creates minimal valid BAM file)
    samtools view -bS header.sam > ${prefix}_all_reads.bam 2>/dev/null || touch ${prefix}_all_reads.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.5
        samtools: 1.24
    END_VERSIONS
    """
}
