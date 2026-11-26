// argNorm: maps each hit of an ARG table to the CARD Antibiotic Resistance
// Ontology (ARO), with its drug and drug-class AROs. Its mapping tables ship
// inside the package (CARD 4.0.0 RGI hits on the source database's sequences),
// so no database is staged. Called as ARGNORM_SARG (ARGs-OAP stage two's
// blastout.filtered.txt) and ARGNORM_DEEPARG (DEEPARG_CONTIGS' .mapping.ARG).
// A local module, not nf-core's argnorm (0.8.0): the ARGs-OAP 4.0.0 table has a
// header row that argNorm, written for 2.3's headerless table, would read as a
// hit, and the input with no hit needs its own output (roadmap Phase 27b).
process ARGNORM {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::argnorm=1.1.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/argnorm:1.1.0--pyhdfd78af_0' :
        'quay.io/biocontainers/argnorm:1.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(hits)
    val tool

    output:
    tuple val(meta), path("*.argnorm.tsv"), emit: tsv
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(tool in ['argsoap', 'deeparg'])) {
        error "ARGNORM: tool must be argsoap or deeparg, not '${tool}'"
    }
    """
    # ARGs-OAP 4.0.0 writes a header row (qseqid sseqid ...); argNorm reads the
    # table without one and takes column 2 (the SARG sequence id), so the header
    # is dropped and the output columns are named 0-17
    if [ "${tool}" = "argsoap" ] && head -n 1 ${hits} | grep -q '^qseqid'; then
        tail -n +2 ${hits} > argnorm_input.tsv
    else
        cp ${hits} argnorm_input.tsv
    fi

    if [ "${tool}" = "argsoap" ] && [ ! -s argnorm_input.tsv ]; then
        # no hit: argNorm cannot parse an empty table, so write its header only
        {
            echo "# argNorm version: \$(argnorm --version)"
            printf '%s\\t' \$(seq 0 17)
            printf 'ARO\\tARO_name\\tCut_Off\\tconfers_resistance_to\\tconfers_resistance_to_names\\tresistance_to_drug_classes\\tresistance_to_drug_classes_names\\n'
        } > ${prefix}.argnorm.tsv
    else
        argnorm ${tool} -i argnorm_input.tsv -o ${prefix}.argnorm.tsv ${args}
    fi
    rm argnorm_input.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        argnorm: \$(argnorm --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "# argNorm version: 1.1.0" > ${prefix}.argnorm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        argnorm: 1.1.0
    END_VERSIONS
    """
}
