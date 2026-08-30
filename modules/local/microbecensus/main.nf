/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MICROBECENSUS Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Average genome size and genome equivalents from host-removed reads
    (design doc Section 4.7, task T5), consumed by AGGREGATE_FUNCTIONS for
    the copies-per-genome-equivalent normalization (Section 6.2).

    The tool is invoked through bin/run_microbe_census_py3fix.py: the pinned
    py3 biocontainer ships an unreleased-upstream-fix gap that crashes the
    stock CLI before any work (upstream issue #32; the py2 image is broken
    differently — Q4 register row has the full story). The wrapper
    monkeypatches the one broken function and delegates; everything else is
    the stock tool. Reads are passed comma-joined (documented paired-end
    form; reads are treated independently, so the singleton file rides
    along). Read length is auto-detected and rounded DOWN to the nearest
    trained model (50..150 by 10, then 175..500); a median below 50 bp makes
    the tool itself exit with an explicit error, per the spec.

    Failure is NON-FATAL by design (Section 4.7): the withName block in
    config/modules.config sets errorStrategy to retry-on-OOM then 'ignore',
    so a failed sample simply contributes no ags.tsv and AGGREGATE_FUNCTIONS
    falls back to TPM-only for it, with an 'unavailable' row in
    ags_and_ge.tsv. The raw key-value output is parsed here into the
    canonical <id>.ags.tsv, rejecting missing/non-numeric/non-positive
    fields and an average genome size outside 0.5-20 Mb (upstream issue #36
    reports occasional silent implausible estimates at zero exit status).

    versions.yml hardcodes 1.1.1: the tool misreports itself as 1.1.0
    (upstream never bumped __version__ for the v1.1.1 release; verified
    against the pinned container) — the container tag is authoritative.

    Input:  tuple val(meta), path(reads)  host-removed reads (R1/R2[/S])
    Output: tuple val(meta), path("*.ags.tsv"), emit: ags
            tuple val(meta), path("*.microbecensus.txt"), emit: raw
            path "versions.yml", emit: versions
----------------------------------------------------------------------------------------
*/

process MICROBECENSUS {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/microbecensus:1.1.1--pyhca03a8a_2' :
        'quay.io/biocontainers/microbecensus:1.1.1--pyhca03a8a_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.ags.tsv")          , emit: ags
    tuple val(meta), path("*.microbecensus.txt"), emit: raw
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = reads instanceof List ? reads.join(',') : "${reads}"
    """
    export TMPDIR=.

    run_microbe_census_py3fix.py \\
        -t ${task.cpus} \\
        ${args} \\
        ${input_reads} \\
        ${prefix}.microbecensus.txt

    # Parse the key-value output ("key:<TAB>value") into the canonical
    # ags.tsv, refusing anything missing, non-numeric, non-positive, or an
    # AGS outside the 0.5-20 Mb plausibility window (header comment)
    awk -F'\\t' -v sample='${prefix}' '
        \$1 == "average_genome_size:" { ags = \$2 }
        \$1 == "genome_equivalents:"  { ge  = \$2 }
        \$1 == "total_bases:"         { tb  = \$2 }
        END {
            num = "^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?\$"
            if (ags !~ num || ge !~ num || tb !~ num) {
                print "MICROBECENSUS: missing or non-numeric field(s) in output (average_genome_size=" ags ", genome_equivalents=" ge ", total_bases=" tb ")" > "/dev/stderr"
                exit 1
            }
            if (ags + 0 <= 0 || ge + 0 <= 0 || tb + 0 <= 0) {
                print "MICROBECENSUS: non-positive estimate(s) (average_genome_size=" ags ", genome_equivalents=" ge ", total_bases=" tb ")" > "/dev/stderr"
                exit 1
            }
            if (ags + 0 < 500000 || ags + 0 > 20000000) {
                print "MICROBECENSUS: average_genome_size " ags " bp outside the 0.5-20 Mb plausibility window; refusing the estimate (upstream issue #36)" > "/dev/stderr"
                exit 1
            }
            print "sample_id\\taverage_genome_size_bp\\tgenome_equivalents\\ttotal_bases"
            print sample "\\t" ags "\\t" ge "\\t" tb
        }' ${prefix}.microbecensus.txt > ${prefix}.ags.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        microbecensus: 1.1.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ags.tsv
    touch ${prefix}.microbecensus.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        microbecensus: 1.1.1
    END_VERSIONS
    """
}
