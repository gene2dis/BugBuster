/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUPERFOCUS Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUPER-FOCUS read-level functional backend (design doc Section 4.6.3, task
    T8b, Q17): translated search of the sample's reads against the SUPER-FOCUS
    DB_90 SEED cluster database (DIAMOND blastx or MMseqs2, selected by
    --superfocus_aligner), then bin/superfocus_function_profile.py sums the
    function-level read counts to SEED subsystem levels 1-3.

    Input handling, each verified on SUPER-FOCUS 1.8 (2026-10-07):
      - only *.fastq/*.fasta/*.fna queries are read (.gz is silently
        dropped) and every query FILE becomes its own count column, so R1, R2
        and the singleton file are decompressed and concatenated into one
        <prefix>.fastq: one column per sample, mates counted as separate
        reads (they sit far apart in the file, so SUPER-FOCUS's grouping of
        consecutive identical read ids never merges a pair)
      - an empty query makes both aligners exit 1 -> a sample without reads
        short-circuits to zero tables without invoking the tool
      - -b points at the database ROOT (containing db/database_PKs.txt and
        db/static/<aligner>/); the 1.8 package ships no db/ folder of its own
    The cluster level is fixed at DB_90 (owner decision, Q17); identity /
    alignment-length / e-value thresholds come from task.ext.args. Alignments
    are deleted (-d); the per-read binning table is not published (one row
    per read hit).

    The aligner is resolved from params.superfocus_aligner in the script
    (FEATURECOUNTS_GENES / RGI_BWT precedent): the composer must know it to
    verify the 'Aligner used:' line of the output.

    versions.yml records the live SUPER-FOCUS / DIAMOND / MMseqs2 versions
    and the database provenance from the DB_VERSION file
    bin/superfocus_db_reformat.sh writes (else 'custom (<aligner> DB_90)').

    Input:  tuple val(meta), path(reads), path(sf_db)
    Output: tuple val(meta), path("*.superfocus_functions.tsv"), emit: functions
            tuple val(meta), path("*.superfocus_summary.tsv"), emit: summary
            tuple val(meta), path("*.superfocus_{all_levels_and_function,subsystem_level_*}.xls"),
                  emit: raw, optional (absent for a sample without reads)
            tuple val(meta), path("*.superfocus.log"), emit: log
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process SUPERFOCUS {
    tag "${meta.id}"
    label 'process_medium'

    // Docker image for every engine: the image's oras:// singularity variant
    // has no gzip (design doc Q20); singularity/apptainer convert the docker
    // image at pull time
    container 'community.wave.seqera.io/library/super-focus_diamond_mmseqs2_unzip_wget:72a2f2b49608ab17'

    input:
    tuple val(meta), path(reads), path(sf_db)

    output:
    tuple val(meta), path("*.superfocus_functions.tsv"), emit: functions
    tuple val(meta), path("*.superfocus_summary.tsv")  , emit: summary
    tuple val(meta), path("*.superfocus_{all_levels_and_function,subsystem_level_*}.xls"), emit: raw, optional: true
    tuple val(meta), path("*.superfocus.log")          , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligner = params.superfocus_aligner
    if ( !(aligner in ['diamond', 'mmseqs2']) ) {
        error "SUPERFOCUS: --superfocus_aligner must be 'diamond' or 'mmseqs2', got '${aligner}'"
    }
    def aligner_db = aligner == 'diamond' ? 'diamond/90_clusters.db.dmnd' : 'mmseqs2/90_clusters.db'
    def read_list = reads instanceof List ? reads.join(' ') : "${reads}"
    """
    if [ ! -s ${sf_db}/db/database_PKs.txt ]; then
        if [ -s ${sf_db}/database_PKs.txt ]; then
            echo "ERROR: ${sf_db} looks like the db/ folder itself; point --custom_superfocus_db at its PARENT (the directory containing db/)" >&2
        else
            echo "ERROR: ${sf_db}/db/database_PKs.txt not found - not a SUPER-FOCUS database root" >&2
        fi
        exit 1
    fi
    if [ ! -s ${sf_db}/db/static/${aligner_db} ]; then
        echo "ERROR: ${sf_db}/db/static/${aligner_db} not found - the database has no ${aligner} DB_90 files (--superfocus_aligner ${aligner})" >&2
        exit 1
    fi

    gzip -cdf ${read_list} > ${prefix}.fastq
    n_lines=\$(wc -l < ${prefix}.fastq)
    if [ \$((n_lines % 4)) -ne 0 ]; then
        echo "ERROR: ${prefix}.fastq has \${n_lines} lines, not a multiple of 4 (expected 4-line FASTQ records)" >&2
        exit 1
    fi
    n_reads=\$((n_lines / 4))

    if [ "\${n_reads}" -eq 0 ]; then
        # an empty query makes diamond/mmseqs exit 1 (verified on 1.8)
        echo "No reads for ${meta.id}: SUPER-FOCUS not run" > ${prefix}.superfocus.log
        table_arg=""
    else
        superfocus \\
            -q ${prefix}.fastq \\
            -dir superfocus_out \\
            -a ${aligner} \\
            -db DB_90 \\
            -b ${sf_db} \\
            -t ${task.cpus} \\
            -tmp ./superfocus_tmp \\
            -n 1 \\
            -d \\
            -l ${prefix}.superfocus.log \\
            ${args} \\
            || { cat ${prefix}.superfocus.log >&2; exit 1; }

        mv superfocus_out/output_all_levels_and_function.xls ${prefix}.superfocus_all_levels_and_function.xls
        for level in 1 2 3; do
            mv superfocus_out/output_subsystem_level_\${level}.xls ${prefix}.superfocus_subsystem_level_\${level}.xls
        done
        table_arg="--table ${prefix}.superfocus_all_levels_and_function.xls"
    fi

    superfocus_function_profile.py \\
        \${table_arg} \\
        --query ${prefix}.fastq \\
        --input-reads \${n_reads} \\
        --aligner ${aligner} \\
        --database 90 \\
        --sample-id ${meta.id} \\
        --prefix ${prefix}

    rm -f ${prefix}.fastq

    if [ -s ${sf_db}/DB_VERSION ]; then
        sf_db_version=\$(head -1 ${sf_db}/DB_VERSION)
    else
        sf_db_version="custom (${aligner} DB_90)"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        superfocus: \$(superfocus -v 2>&1 | sed 's/.*version //')
        diamond: \$(diamond version 2>&1 | sed 's/diamond version //')
        mmseqs2: \$(mmseqs version 2>&1)
        superfocus_aligner: ${aligner}
        superfocus_db: \${sf_db_version}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.superfocus_functions.tsv
    touch ${prefix}.superfocus_summary.tsv
    touch ${prefix}.superfocus_all_levels_and_function.xls
    touch ${prefix}.superfocus_subsystem_level_{1,2,3}.xls
    touch ${prefix}.superfocus.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        superfocus: 1.8
        diamond: 2.2.1
        mmseqs2: 18.8cc5c
        superfocus_aligner: ${params.superfocus_aligner}
        superfocus_db: stub
    END_VERSIONS
    """
}
