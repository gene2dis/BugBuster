/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GALAH_CLUSTER Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run-level dereplication of the refined MAGs (--dereplicate_mags, roadmap
    Phase 18b, B5): galah clusters the refined bins of every assembly of the
    run at --dereplicate_ani (skani) and picks each cluster's representative
    by completeness - 5 x contamination from Binette's CheckM2 reports, the
    choice dRep makes (docs/benchmarks/multisample_binning_2026.md). The
    clusters cover this run's samples only.

    Each MAG is named <assembly>_<bin> (the Bakta file naming); the cluster
    table keeps the assembly and the bin in their own columns, from the
    inputs, since an assembly id may contain '_'.

    Input:
        tuple val(ids), path(bins), path(reports)
            ids:     assembly ids, in the order of bins and reports
            bins:    each assembly's refined bin directory (<id>_binette_bins, *.fa)
            reports: each assembly's CheckM2 report of its refined bins (Name = bin)

    Output:
        clusters:        path("mag_clusters.tsv") - assembly, bin, genome, representative, ani_pct
        representatives: path("representatives") - <assembly>_<bin>.fa of each cluster's representative
        versions:        path(versions.yml)
----------------------------------------------------------------------------------------
*/

process GALAH_CLUSTER {
    tag "dereplicate"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/galah:0.6.2--hfa8f182_0' :
        'quay.io/biocontainers/galah:0.6.2--hfa8f182_0' }"

    input:
    tuple val(ids), path(bins, stageAs: 'bins/*'), path(reports, stageAs: 'reports/*')

    output:
    path "mag_clusters.tsv", emit: clusters
    path "representatives" , emit: representatives
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // the ANI galah clustered at, a column of the cluster table (ext.ani: --dereplicate_ani)
    def ani = task.ext.ani ?: ''
    def bin_list = bins instanceof List ? bins : [bins]
    def report_list = reports instanceof List ? reports : [reports]
    if (ids.size() != bin_list.size() || ids.size() != report_list.size()) {
        error "GALAH_CLUSTER: ${ids.size()} assemblies, ${bin_list.size()} bin directories and ${report_list.size()} reports"
    }
    def add_assemblies = [ids, bin_list, report_list].transpose().collect { id, dir, report ->
        "add_assembly '${id}' '${dir}' '${report}'"
    }.join('\n    ')
    """
    set -euo pipefail

    # genomes/<assembly>_<bin>.fa, one CheckM2 report under those names, and
    # the name -> assembly, bin manifest
    mkdir genomes
    : > genomes.txt
    printf 'genome\\tassembly\\tbin\\n' > manifest.tsv
    rm -f checkm2.tsv
    add_assembly() {
        local id=\$1 dir=\$2 report=\$3 f bin
        [ -e checkm2.tsv ] || head -n 1 "\$report" > checkm2.tsv
        for f in "\$dir"/*.fa; do
            [ -e "\$f" ] || continue
            bin=\$(basename "\$f" .fa)
            if [ -e "genomes/\${id}_\${bin}.fa" ]; then
                echo "ERROR: two MAGs are named \${id}_\${bin}" >&2
                exit 1
            fi
            ln -s "../\$f" "genomes/\${id}_\${bin}.fa"
            echo "genomes/\${id}_\${bin}.fa" >> genomes.txt
            printf '%s\\t%s\\t%s\\n' "\${id}_\${bin}" "\$id" "\$bin" >> manifest.tsv
            awk -F '\\t' -v OFS='\\t' -v b="\$bin" -v g="\${id}_\${bin}" 'NR > 1 && \$1 == b { \$1 = g; print }' "\$report" >> checkm2.tsv
        done
    }
    ${add_assemblies}

    galah cluster \\
        --genome-fasta-list genomes.txt \\
        --checkm2-quality-report checkm2.tsv \\
        --quality-formula completeness-5contamination \\
        ${args} \\
        --threads ${task.cpus} \\
        --output-cluster-definition clusters.tsv

    # assembly, bin, genome, representative, ani_pct (galah writes
    # representative, member paths); a MAG galah's quality filter left out is
    # in no cluster
    mkdir representatives
    awk -F '\\t' -v OFS='\\t' -v ani='${ani}' '
        NR == FNR { if (FNR > 1) { asm[\$1] = \$2; bin[\$1] = \$3 } next }
        function name(p) { sub(/.*\\//, "", p); sub(/[.]fa\$/, "", p); return p }
        { g = name(\$2); r = name(\$1); print asm[g], bin[g], g, r, ani }
    ' manifest.tsv clusters.tsv | sort -t "\$(printf '\\t')" -k4,4 -k3,3 > body.tsv
    { printf 'assembly\\tbin\\tgenome\\trepresentative\\tani_pct\\n'; cat body.tsv; } > mag_clusters.tsv
    cut -f4 body.tsv | sort -u | while read -r r; do
        cp -L "genomes/\${r}.fa" representatives/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: \$(galah --version 2>&1 | sed 's/galah //')
        skani: \$(skani --version 2>&1 | sed 's/skani //')
    END_VERSIONS
    """

    stub:
    """
    mkdir representatives
    printf '>contig_1\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n' > representatives/${ids[0]}_binette_bin1.fa
    printf 'assembly\\tbin\\tgenome\\trepresentative\\tani_pct\\n' > mag_clusters.tsv
    for id in ${ids.join(' ')}; do
        printf '%s\\tbinette_bin1\\t%s_binette_bin1\\t${ids[0]}_binette_bin1\\t${task.ext.ani ?: ''}\\n' "\$id" "\$id" >> mag_clusters.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        galah: 0.6.2
        skani: 0.4.0
    END_VERSIONS
    """
}
