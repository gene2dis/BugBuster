/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN_DBCAN Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CAZy family annotation of the predicted proteins with run_dbcan v5
    (design doc Section 4.3, task T6), protein-mode CAZyme_annotation only:
    DIAMOND vs CAZy, pyHMMER vs dbCAN-HMM, pyHMMER vs dbCAN-sub. The
    consolidated overview.tsv retains the per-tool columns; consensus policy
    is applied downstream in AGGREGATE_FUNCTIONS (params.dbcan_consensus),
    never here. dbCANsub_hmm_results.tsv carries the protein-level substrate
    predictions (dbCAN-sub subfamily mapping); v5 moved CGC-level substrate
    prediction behind CGC identification, which is out of scope.

    The FAA is decompressed with gzip -cdf (Pyrodigal emits .faa.gz; plain
    FASTA passes through unchanged — DeepARG/eggNOG precedent). Empty gene
    sets (header-only FAA from empty-contig samples) short-circuit before
    the tool runs: run_dbcan has no zero-sequence path, so the module emits
    a header-only overview.tsv (exact v5 8-column header incl. the Substrate
    column real output appends, so the AGGREGATE_FUNCTIONS layout guard
    still passes) plus empty companion files (FEATURECOUNTS_GENES
    precedent).

    Output files are renamed to per-sample prefixes: AGGREGATE_FUNCTIONS
    stages every sample's overview into one directory, so the tool's fixed
    output names would collide.

    versions.yml records the live `run_dbcan version` string and the
    database release from the DB_VERSION file that bin/dbcan_db_reformat.sh
    writes into the db dir ('custom' when absent, i.e. --custom_dbcan_db).

    Input:  tuple val(meta), path(proteins), path(dbcan_db)
    Output: tuple val(meta), path("*.overview.tsv"), emit: overview
            tuple val(meta), path("*.dbCANsub_hmm_results.tsv"), emit: substrates
            tuple val(meta), path("*.dbCAN_hmm_results.tsv"), emit: hmm_raw
            tuple val(meta), path("*.diamond.out"), emit: diamond_raw
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process RUN_DBCAN {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.9--pyhdfd78af_0' :
        'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(proteins), path(dbcan_db)

    output:
    tuple val(meta), path("*.overview.tsv")           , emit: overview
    tuple val(meta), path("*.dbCANsub_hmm_results.tsv"), emit: substrates
    tuple val(meta), path("*.dbCAN_hmm_results.tsv")  , emit: hmm_raw
    tuple val(meta), path("*.diamond.out")            , emit: diamond_raw
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${proteins} > ${prefix}_proteins_input.faa

    n_proteins=\$(grep -c '^>' ${prefix}_proteins_input.faa || true)

    if [ "\${n_proteins}" -eq 0 ]; then
        # Empty gene set (header-only FAA from an empty-contig sample):
        # run_dbcan has no zero-sequence path — emit the exact v5 8-column
        # header (verified on real 5.2.9 output, incl. Substrate) so
        # downstream parsing still validates (header comment)
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' \\
            > ${prefix}.overview.tsv
        touch ${prefix}.dbCANsub_hmm_results.tsv
        touch ${prefix}.dbCAN_hmm_results.tsv
        touch ${prefix}.diamond.out
    else
        run_dbcan CAZyme_annotation \\
            --mode protein \\
            --input_raw_data ${prefix}_proteins_input.faa \\
            --db_dir ${dbcan_db} \\
            --output_dir ${prefix}_dbcan_out \\
            --threads ${task.cpus} \\
            ${args}

        mv ${prefix}_dbcan_out/overview.tsv               ${prefix}.overview.tsv
        mv ${prefix}_dbcan_out/dbCANsub_hmm_results.tsv   ${prefix}.dbCANsub_hmm_results.tsv
        mv ${prefix}_dbcan_out/dbCAN_hmm_results.tsv      ${prefix}.dbCAN_hmm_results.tsv
        mv ${prefix}_dbcan_out/diamond.out                ${prefix}.diamond.out
    fi

    rm -f ${prefix}_proteins_input.faa

    if [ -s ${dbcan_db}/DB_VERSION ]; then
        dbcan_db_version=\$(cat ${dbcan_db}/DB_VERSION)
    else
        dbcan_db_version="custom"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        run_dbcan: \$(run_dbcan version | sed 's/dbCAN version: //')
        dbcan_db: \${dbcan_db_version}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.overview.tsv
    touch ${prefix}.dbCANsub_hmm_results.tsv
    touch ${prefix}.dbCAN_hmm_results.tsv
    touch ${prefix}.diamond.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        run_dbcan: 5.2.9
        dbcan_db: stub
    END_VERSIONS
    """
}
