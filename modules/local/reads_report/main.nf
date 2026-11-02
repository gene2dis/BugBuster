process READS_REPORT {

    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88 conda-forge::typst=0.15.1 conda-forge::jinja2=3.1.6 conda-forge::altair=6.3.0 conda-forge::vl-convert-python=1.9.0.post1 conda-forge::font-inter=4.1 conda-forge::pyyaml=6.0.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:2ce015cec686e756' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:6bff5ca407ba1454' }"

    tag "reads_report"
    label 'process_single'

    input:
        path(reports)
        val(args)

    output:
        path("Reads_report.csv"), emit: report
        path("Box_plot_reads.{pdf,svg,png,csv}"), emit: plot, optional: true  // with its plotted values
        path "versions.yml", emit: versions

    script:
        """
        export MPLCONFIGDIR="\$PWD/.mplconfig" XDG_CACHE_HOME="\$PWD/.cache"
        report_unify.py ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python3 --version 2>&1 | sed 's/Python //')
            pandas: \$(python3 -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo unknown)
        END_VERSIONS
	"""

    stub:
        """
        touch Reads_report.csv
        touch Box_plot_reads.pdf Box_plot_reads.svg Box_plot_reads.png Box_plot_reads.csv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: 3.13.18
            pandas: 3.1.0
        END_VERSIONS
        """
}
