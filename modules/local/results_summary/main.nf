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
    path(rgi_alleles, stageAs: 'in_rgi/*')          // <sample>_rgi_bwt.allele_mapping_data.txt; [] unless RGI ran
    path(rgi_kmer_genes, stageAs: 'in_rgi_kmer/*')  // <sample>_rgi_kmer_<k>mer_analysis.gene.txt; [] without RGI k-mer
    path(software_versions, stageAs: 'in_versions/*') // pipeline_info/software_versions.yml
    path(sample_metadata)                           // optional user metadata TSV/CSV; []
    path(card_prevalence)                           // optional WildCARD card_prevalence.txt; []
    path(arg_priority)                              // ARG tier rule table
    path(watchlist)                                 // watchlist taxa
    path(stage_reports, stageAs: 'in_reports/*')   // run-level tables of other stages, by file name (see the case below); []

    output:
    path("results_summary.html"), emit: html
    path("results_summary.xlsx"), emit: xlsx, optional: true
    path("tables/*.tsv")        , emit: tables
    path("figures/*")           , emit: figures, optional: true
    path("README.md")           , emit: readme
    path("README.docx")         , emit: readme_docx
    path("LLM_GUIDE.md")        , emit: llm_guide
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def kmarvel_db_arg  = kmarvel_db ? '--kmarvel-db kmarvel_db' : ''
    def metadata_arg    = sample_metadata ? "--metadata ${sample_metadata}" : ''
    def prevalence_arg  = card_prevalence ? "--card-prevalence ${card_prevalence}" : ''
    // results_summary.py reads a BugBuster results directory; lay the staged
    // files out as one. The 'batch' column is ext.batch_label
    // (--summary_batch_label, or the output directory's name)
    def run_dir = "run/results"
    def batch_label = task.ext.batch_label ?: 'run'
    // the README, LLM_GUIDE.md and the report name the run's output directory, not this staged copy
    def results_path = task.ext.results_path ? file(task.ext.results_path).toAbsolutePath() : run_dir
    """
    mkdir -p ${run_dir}/01_quality_control/summary ${run_dir}/01_quality_control/fastp \\
        ${run_dir}/02_taxonomy/tables ${run_dir}/02_taxonomy/sylph \\
        ${run_dir}/05_arg_prediction/read_level/kmarvel ${run_dir}/pipeline_info \\
        ${run_dir}/05_arg_prediction/read_level/rgi ${run_dir}/05_arg_prediction/read_level/rgi_kmer
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
    for f in in_rgi/*_rgi_bwt.allele_mapping_data.txt; do
        [ -e "\$f" ] || continue
        s=\$(basename "\$f" _rgi_bwt.allele_mapping_data.txt)
        mkdir -p ${run_dir}/05_arg_prediction/read_level/rgi/"\$s"
        ln -s "\$PWD/\$f" ${run_dir}/05_arg_prediction/read_level/rgi/"\$s"/
    done
    for f in in_rgi_kmer/*mer_analysis.gene.txt; do
        [ -e "\$f" ] || continue
        s=\$(basename "\$f" | sed -E 's/_rgi_kmer_[0-9]+mer_analysis[.]gene[.]txt\$//')
        mkdir -p ${run_dir}/05_arg_prediction/read_level/rgi_kmer/"\$s"
        ln -s "\$PWD/\$f" ${run_dir}/05_arg_prediction/read_level/rgi_kmer/"\$s"/
    done
    for f in in_versions/*; do
        [ -e "\$f" ] && ln -s "\$PWD/\$f" ${run_dir}/pipeline_info/software_versions.yml
    done
    # run-level tables of the other stages, laid out at their published paths by file name; a new stage
    # adds one case here and one mix() into ch_summary_reports in main.nf
    for f in in_reports/*; do
        [ -e "\$f" ] || continue
        name=\$(basename "\$f")
        case "\$name" in
            KARGA_norm.csv|KARGVA_norm.csv)
                dest=${run_dir}/05_arg_prediction/read_level/summary ;;
            read_function_abundance.tsv|read_annotated_fraction.tsv|read_sample_summary.tsv)
                dest=${run_dir}/07_functional_annotation/summary ;;
            *.ags.tsv)
                dest=${run_dir}/07_functional_annotation/microbecensus/\${name%.ags.tsv} ;;
            *)
                echo "WARNING: RESULTS_SUMMARY does not read \$name; skipped" >&2
                continue ;;
        esac
        mkdir -p "\$dest"
        ln -s "\$PWD/\$f" "\$dest/"
    done

    results_summary.py \\
        --results ${run_dir} \\
        --batch-labels '${batch_label}' \\
        --results-paths '${results_path}' \\
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
    touch README.md README.docx LLM_GUIDE.md

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
        pandas: 2.2.3
        matplotlib: 3.9.2
        openpyxl: 3.1.5
    END_VERSIONS
    """
}
