process FORMAT_SM_DB {

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/sourmash", pattern: '*'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/kraken", pattern: '*'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/bowtie_index", pattern: '*'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/", pattern: '*'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/", pattern: '*'

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
    publishDir "${params.output}/downloaded_db/", pattern: '*'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/CheckM2_database/", pattern: '*.dmnd'

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
    publishDir "${params.output}/downloaded_db/bowtie_index", pattern: 'phiX_index'

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

    container 'quay.io/ffuentessantander/r_reports:1.1'
    publishDir "${params.output}/downloaded_db/gtdbtk_db", pattern: '*'

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
