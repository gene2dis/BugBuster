/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BOWTIE2_SAMTOOLS Alignment Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Align reads to assembled contigs using Bowtie2 and generate sorted BAM files
    
    Supports both per-sample and co-assembly modes
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
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' :
        'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' }"

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
    
    def is_coassembly = (meta.id == 'coassembly')
    
    if (is_coassembly) {
        // Co-assembly mode: separate R1/R2/singleton files from collected reads
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
            echo "ERROR: Could not identify R1 and R2 files for alignment" >&2
            echo "Files provided: ${reads}" >&2
            exit 1
        fi
        
        # Check if contig file is empty (no sequences after filtering)
        if ! grep -q '^>' "${contigs}"; then
            echo "WARNING: Contig file is empty for ${prefix} - creating empty BAM" >&2
            
            # Create minimal empty BAM with header only
            printf '@HD\tVN:1.6\tSO:coordinate\n' > empty.sam
            samtools view -bS empty.sam > ${prefix}_all_reads.bam
            
            # Generate empty contig tracking report
            printf 'sample_id\t${prefix}\n' > ${prefix}_empty_contig_report.txt
            printf 'stage\tBOWTIE2_SAMTOOLS\n' >> ${prefix}_empty_contig_report.txt
            printf 'empty_contig_detected\tTRUE\n' >> ${prefix}_empty_contig_report.txt
            printf 'timestamp\t%s\n' "\$(date -u +"%Y-%m-%d %H:%M:%S UTC")" >> ${prefix}_empty_contig_report.txt
            printf 'action_taken\tCreated empty BAM file\n' >> ${prefix}_empty_contig_report.txt
            
            # Generate versions file
            printf '"${task.process}":\n' > versions.yml
            printf '    bowtie2: %s\n' "\$(bowtie2 --version 2>&1 | head -n 1 | sed 's/.*version //; s/ /.*/')" >> versions.yml
            printf '    samtools: %s\n' "\$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')" >> versions.yml
            
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
            -1 \$R1_list \\
            -2 \$R2_list \\
            2> ${prefix}_paired.log \\
            | samtools sort \\
                ${args2} \\
                -@ ${task.cpus} \\
                -o ${prefix}_paired.bam -
        
        # Handle singleton reads if present
        if [ -n "\$Singleton_list" ]; then
            bowtie2 \\
                ${args} \\
                -x ${prefix}_index \\
                -p ${task.cpus} \\
                -U \$Singleton_list \\
                2> ${prefix}_singleton.log \\
                | samtools sort \\
                    ${args2} \\
                    -@ ${task.cpus} \\
                    -o ${prefix}_singleton.bam -
            
            # Merge paired and singleton BAMs
            samtools merge \\
                -@ ${task.cpus} \\
                ${prefix}_all_reads.bam \\
                ${prefix}_paired.bam \\
                ${prefix}_singleton.bam
        else
            # No singletons - rename paired BAM
            mv ${prefix}_paired.bam ${prefix}_all_reads.bam
        fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/.*version //; s/ .*//')
            samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')
        END_VERSIONS
        """
    } else {
        // Per-sample assembly mode
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
            printf '@HD\tVN:1.6\tSO:coordinate\n' > empty.sam
            samtools view -bS empty.sam > ${prefix}_all_reads.bam
            
            # Generate empty contig tracking report
            printf 'sample_id\t${prefix}\n' > ${prefix}_empty_contig_report.txt
            printf 'stage\tBOWTIE2_SAMTOOLS\n' >> ${prefix}_empty_contig_report.txt
            printf 'empty_contig_detected\tTRUE\n' >> ${prefix}_empty_contig_report.txt
            printf 'timestamp\t%s\n' "\$(date -u +"%Y-%m-%d %H:%M:%S UTC")" >> ${prefix}_empty_contig_report.txt
            printf 'action_taken\tCreated empty BAM file\n' >> ${prefix}_empty_contig_report.txt
            
            # Generate versions file
            printf '"%s":\n' "${task.process}" > versions.yml
            printf '    bowtie2: %s\n' "\$(bowtie2 --version 2>&1 | head -n 1 | sed 's/.*version //; s/ .*//')" >> versions.yml
            printf '    samtools: %s\n' "\$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')" >> versions.yml
            
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
            
            # Merge paired and singleton BAMs
            samtools merge \\
                -@ ${task.cpus} \\
                ${prefix}_all_reads.bam \\
                ${prefix}_paired.bam \\
                ${prefix}_singleton.bam
        else
            # No singletons - rename paired BAM
            mv ${prefix}_paired.bam ${prefix}_all_reads.bam
        fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version 2>&1 | head -n 1 | sed 's/.*version //; s/ .*//')
            samtools: \$(samtools --version 2>&1 | head -n 1 | sed 's/samtools //')
        END_VERSIONS
        """
    }

    stub:
    prefix = meta.id
    """
    # Create a minimal valid BAM file header
    # Note: In real stub testing, you might want to use samtools to create a proper empty BAM
    echo -e "@HD\tVN:1.6\tSO:coordinate" > header.sam
    echo -e "@SQ\tSN:contig_1\tLN:1000" >> header.sam
    echo -e "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.5.1" >> header.sam
    
    # Convert to BAM (creates minimal valid BAM file)
    samtools view -bS header.sam > ${prefix}_all_reads.bam 2>/dev/null || touch ${prefix}_all_reads.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.1
        samtools: 1.17
    END_VERSIONS
    """
}

process BOWTIE2_SAMTOOLS_DEPTH {
    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    label 'process_medium'

    // BAM files are intermediate - don't publish to final output  
    // publishDir "${params.output}/workflow/${meta.id}/bowtie2_samtools", pattern: '*_all_reads.bam'

    input:
        tuple val(meta), path(reads), path(bins)

    // Output channels for this process

    output:
        tuple val(meta), path("*_all_reads.bam"), emit: reads

    script:
        def prefix = "${meta.id}"

        """
        for bin in ${bins}/*.fa; do
            bin_name=`echo \${bin} | sed -E 's/.+?b//g' | sed 's/.fa//g' | sed 's/in/bin/g'`

	    bowtie2-build \${bin} ${prefix}_bins_index

            bowtie2 \\
            -x ${prefix}_bins_index \\
            -p $task.cpus \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
	    -S ${prefix}_\${bin_name}_paired_reads.sam 2> ${prefix}_bowtie_map.log 

	    samtools sort -@ $task.cpus \\
                  -o ${prefix}_\${bin_name}_paired_reads.bam \\
                  ${prefix}_\${bin_name}_paired_reads.sam

	    if [[ ${reads[2]} == null ]]; then
		    mv ${prefix}_\${bin_name}_paired_reads.bam ${prefix}_\${bin_name}_all_reads.bam
            fi

	    if [[ ${reads[2]} != null ]]; then
	            bowtie2 \\
	            -x ${prefix}_bins_index \\
	            -p $task.cpus \\
	            -U ${reads[2]} \\
	            -S ${prefix}_\${bin_name}_singletons.sam 2> ${prefix}_bowtie_singleton_map.log

		    samtools sort -@ $task.cpus \\
                          -o ${prefix}_\${bin_name}_singletons_reads.bam \\
                          ${prefix}_\${bin_name}_singletons.sam

		    samtools merge ${prefix}_\${bin_name}_all_reads.bam ${prefix}_\${bin_name}_paired_reads.bam ${prefix}_\${bin_name}_singletons_reads.bam
	            rm -f ${prefix}_bowtie_singleton_map.log
                    rm -f ${prefix}_\${bin_name}_singletons.sam
            fi
            rm -f ${prefix}_\${bin_name}_paired_reads.sam
            rm -f ${prefix}_bins_index*
            rm -f ${prefix}_bowtie_map.log
        done
	""" 
}
