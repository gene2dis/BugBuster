process TAXONOMY_PHYLOSEQ {
    tag "$profiler-$db_name"
    label 'process_medium'
    
    conda "conda-forge::python=3.11 conda-forge::pandas=2.0 conda-forge::numpy=1.24 conda-forge::matplotlib=3.7 conda-forge::seaborn=0.12 conda-forge::h5py=3.8 bioconda::biom-format=2.1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jupyter/scipy-notebook:python-3.11' :
        'jupyter/scipy-notebook:python-3.11' }"
    
    publishDir "${params.output}/reports/read_level/taxonomy", mode: params.publish_dir_mode, pattern: '*.{tsv,h5,png}'
    
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
    path("plots/*.png"), emit: plots, optional: true
    path("versions.yml"), emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def plot_levels_arg = plot_levels ?: 'Phylum,Family,Genus,Species'
    def top_n_arg = top_n ?: 10
    def format_arg = task.ext.format ?: 'both'
    """
    # Install biom-format if not present
    pip install --quiet biom-format==2.1.14 2>/dev/null || true
    
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
}
