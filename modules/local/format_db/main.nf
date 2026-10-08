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

process FORMAT_EGGNOG_DB {
    tag "format_eggnog_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download_extensive'

    input:
        val(db)

    output:
        path("eggnog_db")

    script:
        """
        ${params.eggnog_ref_db[params.eggnog_db]["fmtscript"]} $db
        """

    stub:
        """
        mkdir eggnog_db
        touch eggnog_db/eggnog.db eggnog_db/eggnog.db.fieldpresence.bin \\
            eggnog_db/eggnog.db.taxids.bin eggnog_db/eggnog.taxa.db \\
            eggnog_db/eggnog.taxa.db.traverse.pkl eggnog_db/eggnog_proteins.dmnd \\
            eggnog_db/go-basic.obo
        """
}

process FORMAT_DBCAN_DB {
    tag "format_dbcan_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download_extensive'

    input:
        val(db)

    output:
        path("dbcan_db")

    script:
        """
        ${params.dbcan_ref_db[params.dbcan_db]["fmtscript"]} $db
        """

    stub:
        """
        mkdir dbcan_db
        touch dbcan_db/CAZy.dmnd dbcan_db/dbCAN.hmm dbcan_db/dbCAN-sub.hmm \\
            dbcan_db/fam-substrate-mapping.tsv dbcan_db/DB_VERSION
        """
}

process FORMAT_BAKTA_DB {
    tag "format_bakta_db"
    // Runs in the SAME pinned bakta container the BAKTA_BAKTA module uses
    // (not the shared wget image): provisioning must run amrfinder_update —
    // the official DB tarball bundles an AMRFinderPlus DB too old for the
    // container's AMRFinderPlus binary (T7 acceptance finding, 2026-10-06) —
    // and using the identical image guarantees the refreshed DB matches the
    // binary that will consume it. wget/tar/xz verified present in the image
    // (busybox tar; -xJf --strip-components extraction verified live).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.12.1--pyhdfd78af_0' :
        'quay.io/biocontainers/bakta:1.12.1--pyhdfd78af_0' }"

    label 'process_download_extensive'

    input:
        val(db)

    output:
        path("bakta_db")

    script:
        """
        ${params.bakta_ref_db[params.bakta_db]["fmtscript"]} $db
        """

    stub:
        """
        mkdir -p bakta_db/amrfinderplus-db
        touch bakta_db/version.json bakta_db/DB_VERSION
        """
}

process FORMAT_COG_DB {
    tag "format_cog_db"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download'

    input:
        val(db)

    output:
        path("*.def.tab")

    script:
        """
        ${params.cog_ref_db[params.cog_db]["fmtscript"]} $db
        """

    stub:
        """
        touch cog-24.def.tab
        """
}

process FORMAT_WOLTKA_DB {
    tag "format_woltka_db"
    // Runs in the pinned woltka container the WOLTKA_CLASSIFY module uses
    // (not the shared wget image): the md5 files published with WoLr2 cover
    // the UNCOMPRESSED content, so verification needs xz, which the wget
    // image lacks (T7 finding). wget here is busybox (no --tries), hence the
    // retry loop inside bin/woltka_db_reformat.sh
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/woltka:0.1.7--pyhdfd78af_0' :
        'quay.io/biocontainers/woltka:0.1.7--pyhdfd78af_0' }"

    label 'process_download_extensive'

    input:
        val(db)

    output:
        path("woltka_db")

    script:
        """
        ${params.woltka_ref_db[params.woltka_db]["fmtscript"]} $db
        """

    stub:
        """
        mkdir -p woltka_db/databases/bowtie2 woltka_db/proteins \
            woltka_db/function/kegg woltka_db/function/metacyc woltka_db/function/pfam
        touch woltka_db/databases/bowtie2/WoLr2.rev.1.bt2l woltka_db/proteins/coords.txt.xz \
            woltka_db/DB_VERSION
        """
}

process FORMAT_HUMANN_DB {
    tag "format_humann_db"
    // Shared wget image: wget, GNU tar (gzip auto-detected) and md5sum are all
    // the script needs (verified 2026-10-08; the MetaPhlAn tars are plain .tar
    // and only their .pkl / .bt2l members are used, so no bzip2 is required).
    // Downloads ~71 GB with the full ChocoPhlAn (design doc Q2, T8c).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget:1.21.4--5d7af37cfa52d45f' :
        'community.wave.seqera.io/library/wget:1.21.4--c8b4f4320c34b13d' }"

    label 'process_download_extensive'

    input:
        val(key)

    output:
        path("humann_db")

    script:
        def entry = params.humann_ref_db[key]
        """
        ${entry["fmtscript"]} \
            '${key}' \
            '${entry["chocophlan_url"]}' \
            '${entry["uniref_url"]}' \
            '${entry["utility_url"]}' \
            '${entry["metaphlan_index"]}' \
            '${entry["metaphlan_url"]}' \
            '${entry["metaphlan_md5_url"]}' \
            '${entry["metaphlan_bt2_url"]}' \
            '${entry["metaphlan_bt2_md5_url"]}' \
            '${entry["dbversion"]}'
        """

    stub:
        """
        mkdir -p humann_db/chocophlan humann_db/uniref humann_db/utility_mapping humann_db/metaphlan
        touch humann_db/DB_VERSION
        """
}

process FORMAT_SUPERFOCUS_DB {
    tag "format_superfocus_db"
    // Runs in the pinned SUPER-FOCUS image the SUPERFOCUS module uses: it
    // carries GNU wget and unzip (the shared wget image has no unzip), and
    // the archives are prebuilt for the aligners in that same image (design
    // doc Q17). Downloads only the --superfocus_aligner archive.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'oras://community.wave.seqera.io/library/super-focus_diamond_mmseqs2_unzip_wget:fa641269f461e1f8' :
        'community.wave.seqera.io/library/super-focus_diamond_mmseqs2_unzip_wget:72a2f2b49608ab17' }"

    label 'process_download'

    input:
        val(aligner)

    output:
        path("superfocus_db")

    script:
        def entry = params.superfocus_ref_db[params.superfocus_db]
        """
        ${entry["fmtscript"]} \
            ${aligner} \
            '${entry[aligner]["url"]}' \
            '${entry[aligner]["md5"]}' \
            '${entry["pks_url"]}' \
            '${entry["pks_md5"]}' \
            '${entry["dbversion"]}'
        """

    stub:
        def static_file = aligner == 'diamond' ? 'diamond/90_clusters.db.dmnd' : 'mmseqs2/90_clusters.db'
        """
        mkdir -p superfocus_db/db/static/${aligner}
        touch superfocus_db/db/database_PKs.txt superfocus_db/db/static/${static_file} \
            superfocus_db/DB_VERSION
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
