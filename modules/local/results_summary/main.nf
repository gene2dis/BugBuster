process RESULTS_SUMMARY {
    tag "results_summary"
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::pandas=2.2.3 conda-forge::numpy=1.26.4 conda-forge::matplotlib-base=3.9.2 conda-forge::openpyxl=3.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_openpyxl:9bdbbc7916528cdc' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_openpyxl:26abe1de4b35cb08' }"

    input:
    path(reads_report, stageAs: 'in_qc/*')          // Reads_report.csv (TAXONOMY_REPORT's, with the classified fraction, when a profiler ran)
    path(fastp_jsons, stageAs: 'in_fastp/*')        // <sample>.fastp.json; [] without QC
    path(taxonomy_tables, stageAs: 'in_tax/*')      // <label>_{otu,tax}_table.tsv (+ _read_counts.tsv); [] without a profiler
    path(sylph_profiles, stageAs: 'in_sylph/*')     // <sample>_sylph_<db>.tsv; [] unless sylph
    path(kmarvel_tsvs, stageAs: 'in_kmarvel/*')     // <sample>_kmarvel.tsv; [] unless K-MARVEL ran
    path(kmarvel_db, stageAs: 'kmarvel_db')         // K-MARVEL CARD dir (snps.txt); []
    path(software_versions, stageAs: 'in_versions/*') // pipeline_info/software_versions.yml
    path(sample_metadata)                           // optional user metadata TSV/CSV; []
    path(card_prevalence)                           // optional WildCARD card_prevalence.txt; []
    path(arg_priority)                              // ARG tier rule table
    path(watchlist)                                 // watchlist taxa

    output:
    path("results_summary.html"), emit: html
    path("results_summary.xlsx"), emit: xlsx, optional: true
    path("tables/*.tsv")        , emit: tables
    path("figures/*")           , emit: figures, optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def kmarvel_db_arg  = kmarvel_db ? '--kmarvel-db kmarvel_db' : ''
    def metadata_arg    = sample_metadata ? "--metadata ${sample_metadata}" : ''
    def prevalence_arg  = card_prevalence ? "--card-prevalence ${card_prevalence}" : ''
    // results_summary.py reads a BugBuster results directory; lay the staged
    // files out as one. Its name is the 'batch' column (ext.batch_label: the
    // output directory's name), stable across -resume
    def run_dir = "${task.ext.batch_label ?: 'run'}/results"
    """
    mkdir -p ${run_dir}/01_quality_control/summary ${run_dir}/01_quality_control/fastp \\
        ${run_dir}/02_taxonomy/tables ${run_dir}/02_taxonomy/sylph \\
        ${run_dir}/05_arg_prediction/read_level/kmarvel ${run_dir}/pipeline_info
    for f in in_qc/*.csv; do
        [ -e "\$f" ] && ln -s "\$PWD/\$f" ${run_dir}/01_quality_control/summary/Reads_report.csv
    done
    for f in in_fastp/*.fastp.json; do
        [ -e "\$f" ] || continue
        s=\$(basename "\$f" .fastp.json)
        mkdir -p ${run_dir}/01_quality_control/fastp/"\$s"
        ln -s "\$PWD/\$f" ${run_dir}/01_quality_control/fastp/"\$s"/
    done
    for f in in_tax/*.tsv; do
        [ -e "\$f" ] && ln -s "\$PWD/\$f" ${run_dir}/02_taxonomy/tables/
    done
    for f in in_sylph/*.tsv; do
        [ -e "\$f" ] && ln -s "\$PWD/\$f" ${run_dir}/02_taxonomy/sylph/
    done
    for f in in_kmarvel/*_kmarvel.tsv; do
        [ -e "\$f" ] || continue
        s=\$(basename "\$f" _kmarvel.tsv)
        mkdir -p ${run_dir}/05_arg_prediction/read_level/kmarvel/"\$s"
        ln -s "\$PWD/\$f" ${run_dir}/05_arg_prediction/read_level/kmarvel/"\$s"/
    done
    for f in in_versions/*; do
        [ -e "\$f" ] && ln -s "\$PWD/\$f" ${run_dir}/pipeline_info/software_versions.yml
    done

    results_summary.py \\
        --results ${run_dir} \\
        --outdir . \\
        --arg-priority ${arg_priority} \\
        --watchlist ${watchlist} \\
        ${kmarvel_db_arg} \\
        ${metadata_arg} \\
        ${prevalence_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        openpyxl: \$(python -c "import openpyxl; print(openpyxl.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p tables figures
    touch results_summary.html results_summary.xlsx tables/samples_qc.tsv figures/qc_read_fate.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        pandas: 2.2.3
        matplotlib: 3.9.2
        openpyxl: 3.1.5
    END_VERSIONS
    """
}
