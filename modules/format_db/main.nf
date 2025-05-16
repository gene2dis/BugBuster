process FORMAT_SM_DB {

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_download'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_download'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_download'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_download'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'

    label 'process_download'

    input:
        val(db)

    output:
        path("*")

    script:
        """
        ${params.taxonomy_files[params.taxdump_files]["fmtscript"]} $db
        """
}
