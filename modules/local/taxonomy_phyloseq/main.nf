process TAXONOMY_PHYLOSEQ {
    tag "$profiler-$db_name"
    label 'process_medium'
    
    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:cd89ac4e30745df7' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:a1f6f67b03eef1d9' }"
    
    input:
    path(input_files)
    val(profiler)
    val(db_name)
    val(plot_levels)
    val(top_n)
    
    output:
    path("*_otu_table.tsv"), emit: otu_table
    path("*_tax_table.tsv"), emit: tax_table
    path("*_sample_metadata.tsv"), emit: sample_metadata
    path("*_phyloseq_data.h5"), emit: phyloseq_h5, optional: true
    path("*_read_counts.tsv"), emit: read_counts, optional: true  // sylph: estimated read pairs per genome
    path("*.png"), emit: plots, optional: true
    path("versions.yml"), emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def plot_levels_arg = plot_levels ?: 'Phylum,Family,Genus,Species'
    def top_n_arg = top_n ?: 10
    def format_arg = task.ext.format ?: 'both'
    """
    taxonomy_phyloseq.py \\
        --profiler ${profiler} \\
        --input-files ${input_files} \\
        --db-name ${db_name} \\
        --output-dir . \\
        --format ${format_arg} \\
        --plot-levels ${plot_levels_arg} \\
        --top-n ${top_n_arg} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        h5py: \$(python -c "import h5py; print(h5py.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${profiler}_${db_name}_otu_table.tsv
    touch ${profiler}_${db_name}_tax_table.tsv
    touch ${profiler}_${db_name}_sample_metadata.tsv
    touch ${profiler}_${db_name}_phyloseq_data.h5
    if [ "${profiler}" = "sylph" ]; then
        touch ${profiler}_${db_name}_read_counts.tsv
    fi
    touch versions.yml
    """
}
