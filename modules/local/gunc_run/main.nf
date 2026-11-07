/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GUNC_RUN Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Chimerism check of an assembly's refined bins with GUNC (roadmap Phase 16,
    G1): Prodigal genes of each bin mapped with DIAMOND to GUNC's reference,
    and the clade separation score at the taxonomic level where the genes
    split the most (the maxCSS level). It finds chimeras across distant
    lineages that CheckM2's completeness and contamination can miss. A bin is
    one *.fa file; its genome name in the table is the file name without .fa,
    the same as CheckM2's Name (GUNC itself cuts a name at its first '.fa',
    so the bins go in as bin_<n>.fa and the names are mapped back).

    Input:
        tuple val(meta), path(bins)  - refined bin directory (*.fa)
        path(db)                     - GUNC DIAMOND reference (.dmnd, under
                                       GUNC's file name: GUNC reads which
                                       reference it is from the name)

    Output:
        maxcss:   tuple val(meta), path("<id>_gunc_maxCSS_level.tsv")
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process GUNC_RUN {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gunc:1.1.1--pyhdfd78af_0' :
        'quay.io/biocontainers/gunc:1.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bins)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_gunc_maxCSS_level.tsv"), emit: maxcss
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    set -euo pipefail

    # each bin under a neutral name, mapped back after the run: GUNC names a
    # genome by the part of its file name before the first '.fa', which a
    # sample id can contain ('soil.fall21_metabat_bin.1.fa')
    mkdir -p gunc_in gunc_out gunc_tmp
    n=0
    for bin in ${bins}/*.fa; do
        n=\$((n + 1))
        ln -s "\$(readlink -f "\$bin")" gunc_in/bin_\$n.fa
        printf 'bin_%s\t%s\n' "\$n" "\$(basename "\$bin" .fa)" >> bin_names.tsv
    done
    gunc run \\
        --input_dir gunc_in \\
        --file_suffix .fa \\
        --db_file ${db} \\
        --threads ${task.cpus} \\
        --temp_dir gunc_tmp \\
        --out_dir gunc_out \\
        ${args}

    # GUNC.<reference>.maxCSS_level.tsv: the reference GUNC read from the
    # database file name
    table=\$(ls gunc_out/GUNC.*.maxCSS_level.tsv)
    reference=\$(basename "\$table" .maxCSS_level.tsv)
    reference=\${reference#GUNC.}
    awk -F '\\t' -v OFS='\\t' 'NR == FNR { name[\$1] = \$2; next }
        FNR > 1 { if (!(\$1 in name)) { print "GUNC_RUN: unknown genome " \$1 > "/dev/stderr"; exit 1 } \$1 = name[\$1] }
        { print }' bin_names.tsv "\$table" > ${prefix}_gunc_maxCSS_level.tsv
    rm -rf gunc_in gunc_out gunc_tmp bin_names.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$(gunc --version 2>&1)
        diamond: \$(diamond --version 2>&1 | sed 's/^diamond version //')
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/^Prodigal V\\([0-9.]*\\).*/\\1/p')
        gunc_db: \$reference
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    printf 'genome\\tn_genes_called\\tn_genes_mapped\\tn_contigs\\ttaxonomic_level\\tproportion_genes_retained_in_major_clades\\tgenes_retained_index\\tclade_separation_score\\tcontamination_portion\\tn_effective_surplus_clades\\tmean_hit_identity\\treference_representation_score\\tpass.GUNC\\n' > ${prefix}_gunc_maxCSS_level.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: 1.1.1
        diamond: 2.1.24
        prodigal: 2.6.3
        gunc_db: progenomes_2.1
    END_VERSIONS
    """
}
