process PHYLOSEQ_CONVERTER {
    tag "$db_name"
    label 'process_single'
    
    conda "bioconda::bioconductor-phyloseq=1.54.0"
    // phyloseq 1.54.0 biocontainer, pinned for every engine: depot.galaxyproject.org
    // has no 1.54 singularity image (CONTRIBUTING)
    container 'quay.io/biocontainers/bioconductor-phyloseq:1.54.0--r45hdfd78af_0'
    
    input:
    path(otu_table)
    path(tax_table)
    path(sample_metadata)
    val(db_name)
    
    output:
    path("*.RDS"), emit: phyloseq_rds
    path("versions.yml"), emit: versions
    
    when:
    params.create_phyloseq_rds.toString().toBoolean()
    
    script:
    """
    tables_to_phyloseq_simple.R \\
        ${otu_table} \\
        ${tax_table} \\
        ${sample_metadata} \\
        ${db_name}_phyloseq.RDS
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-phyloseq: \$(Rscript -e "cat(as.character(packageVersion('phyloseq')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${db_name}_phyloseq.RDS
    touch versions.yml
    """
}
