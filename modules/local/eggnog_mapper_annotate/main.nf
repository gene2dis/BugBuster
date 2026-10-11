/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EGGNOG_MAPPER_ANNOTATE Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eggNOG-mapper stage 2 of 2: orthology-transfer annotation of the seed
    orthologs from EGGNOG_MAPPER_SEARCH (-m no_search --annotate_hits_table).
    v3 note: the annotations file has 22 columns (Description was dropped) and
    annotation_confidence is a positional code — downstream parsing must be
    version-aware (design doc Section 4.2).

    emapper.py runs through bin/emapper_ogs_fix.py, which patches the beta6
    OG-string parser for eggNOG 7 OG names containing '|' (upstream #620,
    design doc Q23) and stops if the installed method is not the beta6 one.

    Container: eggNOG-mapper v3 is a beta with no bioconda/biocontainer or
    docker image — upstream publishes only an Apptainer .sif, used under
    singularity/apptainer. Every other engine uses the pipeline-built image of
    the same beta6 (containers/eggnog-mapper/Dockerfile, published to GHCR by
    .github/workflows/eggnog-image.yml; pinned by digest — re-pin after a
    rebuild). Re-pin both to the biocontainer when v3.0.0 final ships (design
    doc Q11).

    Input:
        tuple val(meta), path(seed_orthologs), path(eggnog_db)

    Output:
        annotations: tuple val(meta), path("*.emapper.annotations")
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process EGGNOG_MAPPER_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://data.cgmlab.org/eggnog-mapper/emapper-3.0/eggnog-mapper-3.0.0-beta6.sif' :
        'ghcr.io/gene2dis/bugbuster-eggnog-mapper@sha256:4082fbe1ca8be9adcc6653e986cf227a0f725cafdeb4fff9b2b5fb3777f8f6b8' }"

    input:
    tuple val(meta), path(seed_orthologs), path(eggnog_db)

    output:
    tuple val(meta), path("*.emapper.annotations"), emit: annotations
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    emapper_ogs_fix.py \\
        -m no_search \\
        --annotate_hits_table ${seed_orthologs} \\
        ${args} \\
        --cpu ${task.cpus} \\
        --data_dir ${eggnog_db} \\
        --output ${prefix}

    # Keep the full version string — beta suffixes like 3.0.0-beta6 matter for
    # the version-aware parsing downstream (design doc Section 4.2). The
    # database version comes from eggnog.db's own version table, so it is
    # correct for custom databases too.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | grep -o 'emapper-[0-9][^ /]*' | head -1 | sed 's/emapper-//')
        eggnog_db: \$(python3 -c "import sqlite3; print(sqlite3.connect('${eggnog_db}/eggnog.db').execute('SELECT version FROM version').fetchone()[0])" 2>/dev/null || echo unknown)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emapper.annotations

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: 3.0.0-beta6
        eggnog_db: stub
    END_VERSIONS
    """
}
