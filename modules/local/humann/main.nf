/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HUMANN Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HUMAnN 4.0.0a2 read-level functional backend (design doc Section 4.6.1,
    task T8c, Q2): MetaPhlAn 4.1.2 prescreen, nucleotide search against the
    prescreened ChocoPhlAn pangenomes, translated search of the rest against
    UniRef90 (EC-filtered), gene families and MetaCyc pathways in RPK; then
    bin/humann_function_profile.py regroups the gene families to KO / EC and
    takes the pathway abundance, with an exact version/layout guard.

    Pinned stack (Q2, verified 2026-10-08): humann 4.0.0a2 (PyPI) +
    metaphlan 4.1.2 + MetaPhlAn DB mpa_vOct22_CHOCOPhlAnSGB_202403 - the only
    combination 4.0.0a2 accepts - in a Seqera Containers image (no HUMAnN 4
    conda package or biocontainer exists). Docker works.

    Behaviours of 4.0.0a2 handled here, each verified on the pinned image:
      - HUMAnN checks the MetaPhlAn DB tag by running `metaphlan --version`
        without a database option, so the DB is only found through
        METAPHLAN_DB_DIR (exported); --metaphlan-options REPLACES HUMAnN's
        default '-t rel_ab_w_read_stats', so that flag is passed explicitly.
        MetaPhlAn 4.1.2 names its database option --bowtie2db (--db_dir only
        exists from 4.2 and is rejected as an unrecognized argument - found
        on the real-database acceptance run)
      - --utility-database does NOT redirect the MetaCyc pathway files, so
        --pathways-database names both files from the database root
      - an empty input exits 1 ('Unable to determine the input file format')
        -> a sample without reads short-circuits to empty tables without
        invoking the tool
      - R1, R2 and the singleton file are concatenated into one FASTQ under a
        fixed name (mates are separate reads, the HUMAnN convention). Sample
        ids containing 's__' or 't__' are rejected at the input check (a
        profile header line holding both crashes 4.0.0a2's prescreen parser,
        Q2 bug 2)
    The MetaPhlAn profile is written next to the outputs (not into the temp
    dir), so it survives --remove-temp-output.

    versions.yml records the live humann / MetaPhlAn / DIAMOND / Bowtie2
    versions and the database provenance from the DB_VERSION file
    bin/humann_db_reformat.sh writes (else 'custom').

    Input:  tuple val(meta), path(reads), path(humann_db)
    Output: tuple val(meta), path("*.humann_functions.tsv"), emit: functions
            tuple val(meta), path("*.humann_summary.tsv"), emit: summary
            tuple val(meta), path("*_{2_genefamilies,3_reactions,4_pathabundance}.tsv"),
                  emit: raw, optional (absent for a sample without reads)
            tuple val(meta), path("*_1_metaphlan_profile.tsv"), emit: profile,
                  optional (absent without reads or with --bypass-prescreen)
            tuple val(meta), path("*_0.log"), emit: log, optional
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process HUMANN {
    tag "${meta.id}"
    label 'process_high'

    // Docker image for every engine: the image's oras:// singularity variant
    // has no gzip (design doc Q20); singularity/apptainer convert the docker
    // image at pull time
    container 'community.wave.seqera.io/library/python_metaphlan_diamond_bowtie2_pruned:386bc0e5651a4c44'

    input:
    tuple val(meta), path(reads), path(humann_db)

    output:
    tuple val(meta), path("*.humann_functions.tsv")                              , emit: functions
    tuple val(meta), path("*.humann_summary.tsv")                                , emit: summary
    tuple val(meta), path("*_{2_genefamilies,3_reactions,4_pathabundance}.tsv")  , emit: raw, optional: true
    tuple val(meta), path("*_1_metaphlan_profile.tsv")                           , emit: profile, optional: true
    tuple val(meta), path("*_0.log")                                             , emit: log, optional: true
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // the one MetaPhlAn database HUMAnN 4.0.0a2 accepts (design doc Q2)
    def mpa_index = 'mpa_vOct22_CHOCOPhlAnSGB_202403'
    def read_list = reads instanceof List ? reads.join(' ') : "${reads}"
    """
    for sub in chocophlan uniref utility_mapping metaphlan; do
        if [ ! -d ${humann_db}/\${sub} ]; then
            echo "ERROR: ${humann_db}/\${sub}/ not found - not a HUMAnN database root (chocophlan/, uniref/, utility_mapping/, metaphlan/; see --custom_humann_db)" >&2
            exit 1
        fi
    done
    if [ ! -e ${humann_db}/metaphlan/${mpa_index}.pkl ]; then
        echo "ERROR: ${humann_db}/metaphlan/${mpa_index}.pkl not found - HUMAnN 4.0.0a2 accepts only the MetaPhlAn database ${mpa_index}" >&2
        exit 1
    fi

    gzip -cdf ${read_list} > reads.fastq
    n_lines=\$(wc -l < reads.fastq)
    if [ \$((n_lines % 4)) -ne 0 ]; then
        echo "ERROR: the concatenated reads have \${n_lines} lines, not a multiple of 4 (expected 4-line FASTQ records)" >&2
        exit 1
    fi
    n_reads=\$((n_lines / 4))

    if [ "\${n_reads}" -eq 0 ]; then
        # an empty input makes humann exit 1 (verified on 4.0.0a2)
        table_args=""
    else
        export METAPHLAN_DB_DIR=\$(cd ${humann_db}/metaphlan && pwd -P)
        humann \\
            --input reads.fastq \\
            --input-format fastq \\
            --output . \\
            --output-basename ${prefix} \\
            --threads ${task.cpus} \\
            --nucleotide-database ${humann_db}/chocophlan \\
            --protein-database ${humann_db}/uniref \\
            --utility-database ${humann_db}/utility_mapping \\
            --pathways-database ${humann_db}/utility_mapping/metacyc_reactions_level4ec_only.uniref.bz2,${humann_db}/utility_mapping/metacyc_pathways_structured_filtered_v24_subreactions \\
            --count-normalization RPKs \\
            --metaphlan-options "-t rel_ab_w_read_stats --bowtie2db \${METAPHLAN_DB_DIR} --index ${mpa_index} --nproc ${task.cpus}" \\
            --remove-temp-output \\
            ${args}
        table_args="--genefamilies ${prefix}_2_genefamilies.tsv --pathabundance ${prefix}_4_pathabundance.tsv"
    fi

    humann_function_profile.py \\
        \${table_args} \\
        --utility-db ${humann_db}/utility_mapping \\
        --basename ${prefix} \\
        --input-reads \${n_reads} \\
        --sample-id ${meta.id} \\
        --prefix ${prefix}

    rm -f reads.fastq

    if [ -s ${humann_db}/DB_VERSION ]; then
        humann_db_version=\$(head -1 ${humann_db}/DB_VERSION)
    else
        humann_db_version="custom"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: \$(humann --version 2>&1 | sed 's/^humann v//')
        metaphlan: \$(metaphlan --version 2>&1 | head -1 | awk '{print \$3}')
        diamond: \$(diamond version 2>&1 | sed 's/diamond version //')
        bowtie2: \$(bowtie2 --version 2>/dev/null | head -1 | awk '{print \$NF}')
        humann_db: \${humann_db_version}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.humann_functions.tsv
    touch ${prefix}.humann_summary.tsv
    touch ${prefix}_2_genefamilies.tsv ${prefix}_3_reactions.tsv ${prefix}_4_pathabundance.tsv
    touch ${prefix}_1_metaphlan_profile.tsv ${prefix}_0.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        humann: 4.0.0.alpha.2
        metaphlan: 4.1.2
        diamond: 2.1.24
        bowtie2: 2.5.5
        humann_db: stub
    END_VERSIONS
    """
}
