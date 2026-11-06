process RESULTS_SUMMARY {
    tag "results_summary"
    label 'process_single'

    conda "conda-forge::python=3.13 conda-forge::pandas=3.1.0 conda-forge::numpy=2.5.4 conda-forge::matplotlib-base=3.11.2 conda-forge::seaborn-base=0.13.2 conda-forge::h5py=3.16.0 conda-forge::openpyxl=3.1.5 conda-forge::biopython=1.88 conda-forge::typst=0.15.1 conda-forge::jinja2=3.1.6 conda-forge::altair=6.3.0 conda-forge::vl-convert-python=1.9.0.post1 conda-forge::font-inter=4.1 conda-forge::pyyaml=6.0.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:2ce015cec686e756' :
        'community.wave.seqera.io/library/python_pandas_numpy_matplotlib-base_pruned:6bff5ca407ba1454' }"

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
    // assets/report (Typst and Jinja templates, styles, references) and a hash of its files' contents, which
    // makes an edited template rerun the task on -resume (a directory is cached by its own path, size and mtime)
    tuple path(report_assets, stageAs: 'report_assets'), val(report_assets_hash)
    path(stage_reports, stageAs: 'in_reports/*')   // run-level tables of other stages, by file name (see the case below); []

    output:
    path("results_summary.html"), emit: html        // a link to report/index.html
    path("report")              , emit: site        // the multi-page HTML report
    path("results_summary.pdf") , emit: pdf, optional: true  // none when typst fails (a warning) or with --no-pdf
    path("results_summary.xlsx"), emit: xlsx, optional: true
    path("tables/*.tsv")        , emit: tables
    path("figures/*")           , emit: figures, optional: true
    path("README.md")           , emit: readme
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
    // per_sample/<sample>/ or coassembly/ for the bin depth tables (ext.assembly_mode: --assembly_mode)
    def bin_depth_dir = task.ext.assembly_mode == 'assembly' ? 'per_sample' : 'coassembly'
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
    # coassembly_groups.tsv (--coassembly_group_by): sample -> co-assembly group, staged with the stage tables
    groups=in_reports/coassembly_groups.tsv
    is_group() { [ -e "\$groups" ] && cut -f2 "\$groups" | tail -n +2 | grep -qxF -- "\$1"; }
    group_of() { if [ -e "\$groups" ]; then awk -F '\t' -v s="\$1" 'NR > 1 && \$1 == s { print \$2 }' "\$groups"; fi; }
    # an assembly's folder: coassembly, coassembly/<group>, or per_sample/<sample>
    asm_dir() {
        if [ "\$1" = coassembly ]; then echo coassembly
        elif is_group "\$1"; then echo "coassembly/\$1"
        else echo "per_sample/\$1"; fi
    }
    # an assembly's folder under the id-keyed stage folders (deeparg/, mags/): <id>, or coassembly/<group>
    id_dir() { if is_group "\$1"; then echo "coassembly/\$1"; else echo "\$1"; fi; }
    for f in in_reports/*; do
        [ -e "\$f" ] || continue
        name=\$(basename "\$f")
        case "\$name" in
            args_oap_long.tsv|args_oap_cells.tsv)
                dest=${run_dir}/05_arg_prediction/read_level/args_oap_summary ;;
            read_function_abundance.tsv|read_annotated_fraction.tsv|read_sample_summary.tsv)
                dest=${run_dir}/07_functional_annotation/summary ;;
            function_abundance.tsv|annotated_fraction.tsv|gene_abundance.tsv)
                dest=${run_dir}/07_functional_annotation/summary ;;
            *.ags.tsv)
                dest=${run_dir}/07_functional_annotation/microbecensus/\${name%.ags.tsv} ;;
            *.featureCounts.txt.summary)
                dest=${run_dir}/07_functional_annotation/gene_abundance/\${name%.featureCounts.txt.summary} ;;
            Contig_tax_and_arg_prediction.tsv)
                dest=${run_dir}/05_arg_prediction/contig_level/summary ;;
            *_contigs_deep_arg.out.mapping.ARG)
                dest=${run_dir}/05_arg_prediction/contig_level/deeparg/\$(id_dir "\${name%_contigs_deep_arg.out.mapping.ARG}") ;;
            Arg_cluster_*_cluster.tsv)
                dest=${run_dir}/05_arg_prediction/bin_level/clustering ;;
            contig_filtering_summary.txt|coassembly_groups.tsv)
                dest=${run_dir}/pipeline_info ;;
            *_contig.stats)
                id=\${name%_contig.stats}
                dest=${run_dir}/03_assembly/\$(asm_dir "\$id") ;;
            *_bin_depth.tsv)
                # <sample>_bin_depth.tsv: next to the sample's own assembly in per-sample mode, else
                # its co-assembly's
                s=\${name%_bin_depth.tsv}
                g=\$(group_of "\$s")
                if [ -n "\$g" ]; then
                    dest=${run_dir}/04_binning/coassembly/\$g/coverage
                elif [ "${bin_depth_dir}" = per_sample ]; then
                    dest=${run_dir}/04_binning/per_sample/\$s/coverage
                else
                    dest=${run_dir}/04_binning/coassembly/coverage
                fi ;;
            *_gtdbtk_bac120.tsv|*_gtdbtk_ar53.tsv)
                id=\${name%_gtdbtk_*}
                dest=${run_dir}/04_binning/\$(asm_dir "\$id")/taxonomy/gtdbtk ;;
            *_quality_report.tsv)
                # <id>_<binner>_quality_report.tsv; binner names have no underscore
                id=\${name%_quality_report.tsv}
                id=\${id%_*}
                dest=${run_dir}/04_binning/\$(asm_dir "\$id")/quality/checkm2 ;;
            *_binette_bin*.tsv|*_binette_bin*.txt)
                # Bakta's <id>_<bin>.tsv / .txt (it runs on Binette's refined bins only)
                dest=${run_dir}/07_functional_annotation/mags/\$(id_dir "\${name%_binette_bin*}") ;;
            *)
                echo "WARNING: RESULTS_SUMMARY does not read \$name; skipped" >&2
                continue ;;
        esac
        mkdir -p "\$dest"
        ln -s "\$PWD/\$f" "\$dest/"
    done

    # matplotlib's and fontconfig's caches in the task directory (a read-only image, maybe no home)
    export MPLCONFIGDIR="\$PWD/.mplconfig" XDG_CACHE_HOME="\$PWD/.cache"
    results_summary.py \\
        --results ${run_dir} \\
        --batch-labels '${batch_label}' \\
        --outdir . \\
        --arg-priority ${arg_priority} \\
        --watchlist ${watchlist} \\
        --report-assets report_assets \\
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
        typst: \$(typst --version | sed 's/typst //; s/ .*//')
        jinja2: \$(python -c "import jinja2; print(jinja2.__version__)")
        altair: \$(python -c "import altair; print(altair.__version__)")
        vl-convert: \$(python -c "import vl_convert; print(vl_convert.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p tables figures report/assets
    touch results_summary.html results_summary.pdf results_summary.xlsx tables/samples_qc.tsv
    touch figures/qc_read_fate.png figures/qc_read_fate.pdf figures/qc_read_fate.svg figures/qc_read_fate.csv
    touch report/index.html report/methods.html report/assets/site.css README.md LLM_GUIDE.md

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.13.18
        pandas: 3.1.0
        matplotlib: 3.11.2
        openpyxl: 3.1.5
        typst: 0.15.1
        jinja2: 3.1.6
        altair: 6.3.0
        vl-convert: 1.9.0
    END_VERSIONS
    """
}
