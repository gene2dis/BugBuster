process FORMAT_KRAKEN_DB {
    tag "format_kraken_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.kraken_ref_db[params.kraken2_db]["fmtscript"]} $db
        """

    stub:
        """
        mkdir kraken_db
        touch kraken_db/hash.k2d kraken_db/opts.k2d kraken_db/taxo.k2d
        """
}

process FORMAT_NT_BLAST_DB {
    tag "format_blast_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download_extensive'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.blast_ref_db[params.blast_db]["fmtscript"]} $db
        """

    stub:
        """
        touch nt.00.nin nt.00.nhr nt.00.nsq nt.nal
        """
}

process FORMAT_TAXDUMP_FILES {
    tag "format_taxdump"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.taxonomy_files[params.taxdump_files]["fmtscript"]} $db
        """

    stub:
        """
        touch nodes.dmp names.dmp nucl_gb.accession2taxid
        """
}

process DOWNLOAD_DEEPARG_DB {

    container 'quay.io/ffuentessantander/deeparg:1.0.4'

    label 'process_download'

    output:
        path("deeparg_db"), emit: deeparg_db
        path "versions.yml", emit: versions

    script:
        """
        deeparg \\
            download_data \\
            -o ./deeparg_db

        # The DeepARG data download is unversioned server-side; record the tool
        # version and the download date as the only available provenance.
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: \$(pip show deeparg 2>/dev/null | sed -n 's/Version: //p' || echo 1.0.4)
            deeparg_db: unversioned, downloaded \$(date -u +%Y-%m-%d)
        END_VERSIONS
        """

    stub:
        """
        mkdir deeparg_db

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeparg: 1.0.4
            deeparg_db: stub
        END_VERSIONS
        """
}

process FORMAT_CHECKM2_DB {
    tag "format_checkm2_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download'

    input:
        val(db)

    output:
        path("*.dmnd")

    script:
        """
        ${params.checkm2_ref_db[params.checkm2_db]["fmtscript"]} $db
        """

    stub:
        """
        touch uniref100.KO.1.dmnd
        """
}

process DOWNLOAD_GTDBTK_DB {
    tag "download_gtdbtk_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.gtdbtk_ref_db[params.gtdbtk_db]["fmtscript"]} $db
        """

    stub:
        """
        mkdir gtdbtk_db
        touch gtdbtk_db/metadata.txt
        """
}

process SOURMASH_TAX_PREPARE {

    container 'quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0'

    label 'process_single'

    input:
        path(tax_file)

    output:
        path("*.sqldb"), emit: tax_db
        path "versions.yml", emit: versions

    script:
        """
        sourmash tax prepare -t ${tax_file} -o taxonomy.sqldb -F sql

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sourmash: \$(sourmash --version 2>&1 | sed 's/sourmash //')
        END_VERSIONS
        """

    stub:
        """
        touch taxonomy.sqldb

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sourmash: 4.8.11
        END_VERSIONS
        """
}
