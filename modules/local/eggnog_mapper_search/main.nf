/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EGGNOG_MAPPER_SEARCH Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    eggNOG-mapper stage 1 of 2: DIAMOND search of predicted proteins against
    the eggNOG 7 database, annotation deferred (--no_annot). Splitting search
    from annotation makes the expensive step independently retryable
    (design doc Section 4.2).

    Container: eggNOG-mapper v3 is a beta with no bioconda/biocontainer or
    docker image — upstream publishes only an Apptainer .sif, used under
    singularity/apptainer. Every other engine uses the pipeline-built image of
    the same beta6 (containers/eggnog-mapper/Dockerfile, published to GHCR by
    .github/workflows/eggnog-image.yml; pinned by digest — re-pin after a
    rebuild). Re-pin both to the biocontainer when v3.0.0 final ships (design
    doc Q11).

    Input:
        tuple val(meta), path(proteins), path(eggnog_db)
        proteins: predicted proteins FASTA (gzipped Pyrodigal output accepted)
        eggnog_db: eggNOG 7 data dir (eggnog.db + fieldpresence/taxids caches,
                   eggnog.taxa.db + traverse.pkl, eggnog_proteins.dmnd,
                   go-basic.obo — the emapper-3.0 download_eggnog_data.py layout)

    Output:
        seed_orthologs: tuple val(meta), path("*.emapper.seed_orthologs")
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process EGGNOG_MAPPER_SEARCH {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://data.cgmlab.org/eggnog-mapper/emapper-3.0/eggnog-mapper-3.0.0-beta6.sif' :
        'ghcr.io/gene2dis/bugbuster-eggnog-mapper@sha256:4082fbe1ca8be9adcc6653e986cf227a0f725cafdeb4fff9b2b5fb3777f8f6b8' }"

    input:
    tuple val(meta), path(proteins), path(eggnog_db)

    output:
    tuple val(meta), path("*.emapper.seed_orthologs"), emit: seed_orthologs
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Pyrodigal emits gzipped proteins; gzip -cdf also passes plain FASTA
    # through unchanged
    gzip -cdf ${proteins} > ${prefix}_proteins_input.faa

    # emapper exits 1 when the kernel kills diamond for running out of
    # memory, which the memory retry (exit 130-145) does not catch: report
    # such a failure as 137 instead. diamond peaks at ~80 GB against the
    # eggNOG 7 database for a real metagenome (block size 8, one index chunk)
    status=0
    emapper.py \\
        -m diamond \\
        --no_annot \\
        --itype proteins \\
        ${args} \\
        --cpu ${task.cpus} \\
        -i ${prefix}_proteins_input.faa \\
        --data_dir ${eggnog_db} \\
        --output ${prefix} \\
        > emapper.stdout 2> emapper.stderr || status=\$?
    cat emapper.stdout
    cat emapper.stderr >&2
    if [ "\$status" -ne 0 ]; then
        # emapper prints diamond's failure ("Error running diamond: ... Killed") on stdout
        if grep -qE 'Killed|std::bad_alloc|Cannot allocate memory|[Oo]ut of memory|MemoryError' emapper.stdout emapper.stderr; then
            echo "EGGNOG_MAPPER_SEARCH: diamond ran out of memory with ${task.memory}; exiting 137 (retried with more memory unless --max_memory ${params.max_memory} already caps it: then raise --max_memory)" >&2
            exit 137
        fi
        exit "\$status"
    fi

    rm -f ${prefix}_proteins_input.faa emapper.stdout emapper.stderr

    # Keep the full version string — beta suffixes like 3.0.0-beta6 matter for
    # the version-aware parsing downstream (design doc Section 4.2)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$(emapper.py --version 2>&1 | grep -o 'emapper-[0-9][^ /]*' | head -1 | sed 's/emapper-//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.emapper.seed_orthologs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: 3.0.0-beta6
    END_VERSIONS
    """
}
