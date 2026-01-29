process FORMAT_SM_DB {
    tag "format_sourmash_db"
    container 'ubuntu:22.04'

    label 'process_download_single'

    input:
        path(db)

    output:
        path("*")

    script:
        """
        ${params.sourmash_ref_db[params.sourmash_db]["fmtscript"]} 
	"""
}

process FORMAT_KRAKEN_DB {
    tag "format_kraken_db"
    container 'ubuntu:22.04'

    label 'process_download_single'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.kraken_ref_db[params.kraken2_db]["fmtscript"]} $db
        """
}

process FORMAT_BOWTIE_INDEX {
    tag "format_bowtie_index"
    container 'ubuntu:22.04'

    label 'process_download_single'

    input:
        path(db)

    output:
        path("*")

    script:
        """
        ${params.bowtie_ref_host_index[params.host_db]["fmtscript"]}
        """
}

process FORMAT_NT_BLAST_DB {
    tag "format_blast_db"
    container 'ubuntu:22.04'

    label 'process_download_extensive'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.blast_ref_db[params.blast_db]["fmtscript"]} $db
        """
}

process FORMAT_TAXDUMP_FILES {
    tag "format_taxdump"
    container 'ubuntu:22.04'

    label 'process_download_single'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.taxonomy_files[params.taxdump_files]["fmtscript"]} $db
        """
}

process DOWNLOAD_DEEPARG_DB {

    container 'quay.io/ffuentessantander/deeparg:1.0.4'

    label 'process_download_single'

    output:
        path("deeparg_db")

    script:
        """
        deeparg \\
            download_data \\
            -o ./deeparg_db
        """
}

process FORMAT_CHECKM2_DB {
    tag "format_checkm2_db"
    container 'ubuntu:22.04'

    label 'process_download_single'

    input:
        val(db)

    output:
        path("*.dmnd")

    script:
        """
        ${params.checkm2_ref_db[params.checkm2_db]["fmtscript"]} $db
        """
}

process BUILD_PHIX_BOWTIE2_INDEX {

    container 'quay.io/biocontainers/bowtie2:2.5.3--py310ha0a81b8_0'

    label 'process_download_single'

    input:
        path(phiX_fasta)

    output:
        path("phiX_index")

    script:
        """
        mkdir phiX_index
        bowtie2-build ${phiX_fasta} phiX_index/phiX
        """
}

process DOWNLOAD_GTDBTK_DB {
    tag "download_gtdbtk_db"
    container 'ubuntu:22.04'

    label 'process_download_single'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.gtdbtk_ref_db[params.gtdbtk_db]["fmtscript"]} $db
        """
}

process SOURMASH_TAX_PREPARE {

    container 'quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0'

    label 'process_single'

    input:
        path(tax_file)

    output:
        path("*.sqldb"), emit: tax_db

    script:
        """
        sourmash tax prepare -t ${tax_file} -o taxonomy.sqldb -F sql
        """
}
