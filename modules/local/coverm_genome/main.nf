/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COVERM_GENOME Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Abundance of the run's dereplicated MAGs in one sample (--dereplicate_mags,
    roadmap Phase 18b, B5): CoverM genome maps the sample's read pairs to the
    representatives (strobealign; reads at >= 97 % identity, as the binning
    depths) and reports each representative's relative abundance (% of the
    sample's reads, with an 'unmapped' row), mean and trimmed mean depth and
    covered fraction. Every sample of the run against the same set, so the
    values compare across this run's samples.

    Input:
        tuple val(meta), path(reads), path(representatives)
            reads:           R1, R2[, singleton] (the pairs are mapped)
            representatives: directory of <assembly>_<bin>.fa

    Output:
        abundance: tuple val(meta), path("${meta.id}_coverm_genome.tsv") - Genome,
                   Relative Abundance (%), Mean, Trimmed Mean, Covered Fraction
        versions:  path(versions.yml)
----------------------------------------------------------------------------------------
*/

process COVERM_GENOME {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.8.0--h750ce8b_0' :
        'quay.io/biocontainers/coverm:0.8.0--h750ce8b_0' }"

    input:
    tuple val(meta), path(reads), path(representatives)

    output:
    tuple val(meta), path("${meta.id}_coverm_genome.tsv"), emit: abundance
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    if (!reads || reads.size() < 2) {
        error "COVERM_GENOME: paired reads are required for sample ${meta.id}"
    }
    """
    set -euo pipefail

    # CoverM's temporary BAMs in the task directory
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"

    coverm genome \\
        --coupled ${reads[0]} ${reads[1]} \\
        --genome-fasta-directory ${representatives} \\
        --genome-fasta-extension fa \\
        --methods relative_abundance mean trimmed_mean covered_fraction \\
        ${args} \\
        --threads ${task.cpus} \\
        --output-file coverm.tsv

    # One sample per table: the column names without its read file's name
    awk -F '\\t' -v OFS='\\t' 'NR == 1 { print "Genome", "Relative Abundance (%)", "Mean", "Trimmed Mean", "Covered Fraction"; next } { print }' \\
        coverm.tsv > ${prefix}_coverm_genome.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //')
        strobealign: \$(strobealign --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    printf 'Genome\\tRelative Abundance (%%)\\tMean\\tTrimmed Mean\\tCovered Fraction\\nunmapped\\t50\\tNA\\tNA\\tNA\\n' > ${prefix}_coverm_genome.tsv
    for f in ${representatives}/*.fa; do
        printf '%s\\t50\\t10.5\\t10.2\\t0.99\\n' "\$(basename "\$f" .fa)" >> ${prefix}_coverm_genome.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: 0.8.0
        strobealign: 0.17.0
    END_VERSIONS
    """
}
