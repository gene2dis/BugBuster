process PHYLOSEQ_CONVERTER {
    tag "$db_name"
    label 'process_single'
    
    conda "bioconda::bioconductor-phyloseq=1.44.0 conda-forge::r-optparse=1.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.46.0--r43hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-phyloseq:1.46.0--r43hdfd78af_0' }"
    
    input:
    path(otu_table)
    path(tax_table)
    path(sample_metadata)
    val(db_name)
    
    output:
    path("*.RDS"), emit: phyloseq_rds
    path("versions.yml"), emit: versions
    
    when:
    params.create_phyloseq_rds
    
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
}
