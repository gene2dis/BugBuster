/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WOLTKA_CLASSIFY Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Second half of the Woltka read-level functional backend (design doc
    Section 4.6.2, task T8a): `woltka classify` assigns the WOLTKA_ALIGN
    alignments to WoLr2 ORFs by coordinate overlap (>= 80 % of the read
    inside the ORF, Woltka default), then bin/woltka_function_profile.py
    turns the ORF profile into per-function read counts and RPK for the
    ko / ec / cog / pfam / metacyc ontologies. COG ortholog ids are mapped
    to COG functional-category letters with the NCBI COG definitions table
    (cog_def, the contig branch's cog-24.def.tab) and Pfam accessions are
    reported by name, so both ontologies share the contig branch's
    vocabulary (design doc Q14).

    Fixed classify flags, each verified on woltka 0.1.7 (2026-10-07):
      --no-demux    single-file input is otherwise demultiplexed on the
                    first underscore of the READ NAME, silently splitting
                    one sample into bogus samples
      --digits 6    without it Woltka rounds every ORF count to an integer
                    and drops ones that round to 0, discarding the 1/k
                    fractions of multi-hit reads
      --unassigned  report reads left unassigned by --uniq ambiguity
                    (published as a quality signal)
    task.ext.args carries --uniq when params.woltka_uniq is set.

    Functions are composed per ORF in the bin/ script, NOT with
    `woltka collapse`: collapse counts duplicate map entries and every path
    of a chained map, inflating abundances (rationale in the script header).

    An alignment file without records (no read aligned, or an empty input
    sample) short-circuits to a header-only profile: woltka classify raises
    "Alignment file is empty or unreadable" on it (verified live).

    versions.yml records the live woltka version, the database release
    from the DB_VERSION file bin/woltka_db_reformat.sh writes (a custom
    mirror without one records 'custom (<index basename>)') and the COG
    definitions table name.

    Input:  tuple val(meta), path(sam), path(wol_db), path(cog_def)
    Output: tuple val(meta), path("*.woltka_orf.tsv"), emit: orf_profile
            tuple val(meta), path("*.woltka_functions.tsv"), emit: functions
            tuple val(meta), path("*.woltka_summary.tsv"), emit: summary
            tuple val(meta), path("*.woltka_unassigned.tsv"), emit: unassigned
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process WOLTKA_CLASSIFY {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/woltka:0.1.7--pyhdfd78af_0' :
        'quay.io/biocontainers/woltka:0.1.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sam), path(wol_db), path(cog_def)

    output:
    tuple val(meta), path("*.woltka_orf.tsv")       , emit: orf_profile
    tuple val(meta), path("*.woltka_functions.tsv") , emit: functions
    tuple val(meta), path("*.woltka_summary.tsv")   , emit: summary
    tuple val(meta), path("*.woltka_unassigned.tsv"), emit: unassigned
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ "\$(gzip -cdf ${sam} | head -c 1 | wc -c)" -eq 0 ]; then
        # No alignment records: woltka classify cannot parse an empty file
        printf '#FeatureID\\t%s\\n' "${prefix}" > ${prefix}.woltka_orf.tsv
    else
        woltka classify \\
            --input ${sam} \\
            --coords ${wol_db}/proteins/coords.txt.xz \\
            --no-demux \\
            --digits 6 \\
            --unassigned \\
            ${args} \\
            --to-tsv \\
            --output ${prefix}.woltka_orf.tsv
    fi

    woltka_function_profile.py \\
        --profile ${prefix}.woltka_orf.tsv \\
        --db ${wol_db} \\
        --cog-def ${cog_def} \\
        --sample-id ${meta.id} \\
        --prefix ${prefix}

    if [ -s ${wol_db}/DB_VERSION ]; then
        wol_db_version=\$(head -1 ${wol_db}/DB_VERSION)
    else
        index_file=\$(ls ${wol_db}/databases/bowtie2/*.rev.1.bt2l ${wol_db}/databases/bowtie2/*.rev.1.bt2 2>/dev/null | head -1 || true)
        index_name=\$(basename "\${index_file%.rev.1.bt2*}")
        wol_db_version="custom (\${index_name:-unknown})"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        woltka: \$(woltka --version 2>&1 | sed 's/.*version //')
        wol_db: \${wol_db_version}
        cog_def: ${cog_def.name}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.woltka_orf.tsv
    touch ${prefix}.woltka_functions.tsv
    touch ${prefix}.woltka_summary.tsv
    touch ${prefix}.woltka_unassigned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        woltka: 0.1.7
        wol_db: stub
        cog_def: ${cog_def.name}
    END_VERSIONS
    """
}
