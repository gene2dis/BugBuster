/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FEATURECOUNTS_GENES Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-gene read counts (design doc Section 4.4, task T3): featureCounts over
    the shared Pyrodigal gene coordinates, on the reads-vs-contigs BAM. The
    GFF is first converted to SAF by bin/gff2saf.sh, which rewrites gene ids
    from Prodigal's <seqnum>_<genenum> ID attribute to the <contig>_<genenum>
    form used in the FAA headers, so counts join directly against the
    eggNOG-mapper annotations in T4.

    Counting is read-level, not fragment-level (spec 5.2: "Reads assigned by
    featureCounts"): -p is passed WITHOUT --countReadPairs. subread hard-fails
    when -p sees a BAM without paired-flag records (singleton-only sample,
    record-less BAM) and when paired records appear without -p, so the paired
    invocation runs first and falls back to single-end mode on that specific
    error (verified on subread 2.1.1).

    An empty gene set — the empty-contig sample shape, which arrives together
    with a header-only BAM — would make featureCounts fail on a zero-feature
    annotation; it short-circuits to a header-only count table and an all-zero
    summary so counting never crashes (Q1 note in the design doc).

    Multi-mapping policy comes from params.featurecounts_multimap
    (primary | all | none, schema-validated). With the pipeline's Bowtie2
    defaults (one reported alignment per read, no secondary records) the
    settings coincide today; the param keeps the policy explicit as the spec
    requires, and matters if BOWTIE2_SAMTOOLS ext.args ever changes.

    Input:
        tuple val(meta), path(bam), path(gff)
        bam: coordinate-sorted reads-vs-contigs BAM (no index needed)
        gff: Pyrodigal gene coordinates (gzipped accepted)

    Output:
        counts:  tuple val(meta), path("*.featureCounts.txt")
                 (column 6 "Length" is required downstream for TPM — keep it)
        summary: tuple val(meta), path("*.featureCounts.txt.summary")
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process FEATURECOUNTS_GENES {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/subread:2.1.1--h577a1d6_0' :
        'quay.io/biocontainers/subread:2.1.1--h577a1d6_0' }"

    input:
    tuple val(meta), path(bam), path(gff)

    output:
    tuple val(meta), path("*.featureCounts.txt")        , emit: counts
    tuple val(meta), path("*.featureCounts.txt.summary"), emit: summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def multimap_flag = params.featurecounts_multimap == 'primary' ? '--primary' :
                        params.featurecounts_multimap == 'all'     ? '-M' : ''
    """
    gff2saf.sh ${gff} ${prefix}_genes.saf

    if [ \$(tail -n +2 ${prefix}_genes.saf | wc -l) -eq 0 ]; then
        # Empty gene set: emit the header-only table and all-zero summary in
        # featureCounts' own layout so downstream parsing stays uniform
        {
            echo "# Program:featureCounts \$(featureCounts -v 2>&1 | grep featureCounts); empty gene set for ${prefix}, counting skipped"
            printf 'Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t%s\\n' "${bam}"
        } > ${prefix}.featureCounts.txt
        {
            printf 'Status\\t%s\\n' "${bam}"
            for s in Assigned Unassigned_Unmapped Unassigned_Read_Type \\
                     Unassigned_Singleton Unassigned_MappingQuality \\
                     Unassigned_Chimera Unassigned_FragmentLength \\
                     Unassigned_Duplicate Unassigned_MultiMapping \\
                     Unassigned_Secondary Unassigned_NonSplit \\
                     Unassigned_NoFeatures Unassigned_Overlapping_Length \\
                     Unassigned_Ambiguity; do
                printf '%s\\t0\\n' "\$s"
            done
        } > ${prefix}.featureCounts.txt.summary
    else
        # Paired-first with single-end fallback (see header comment)
        set +e
        featureCounts \\
            -F SAF \\
            -a ${prefix}_genes.saf \\
            -T ${task.cpus} \\
            -p ${multimap_flag} \\
            ${args} \\
            -o ${prefix}.featureCounts.txt \\
            ${bam} 2> fc_paired.log
        fc_status=\$?
        set -e
        if [ \$fc_status -ne 0 ]; then
            if grep -q 'No paired-end reads were detected' fc_paired.log; then
                featureCounts \\
                    -F SAF \\
                    -a ${prefix}_genes.saf \\
                    -T ${task.cpus} \\
                    ${multimap_flag} \\
                    ${args} \\
                    -o ${prefix}.featureCounts.txt \\
                    ${bam}
            else
                cat fc_paired.log >&2
                exit \$fc_status
            fi
        else
            cat fc_paired.log >&2
        fi
    fi

    rm -f fc_paired.log ${prefix}_genes.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | sed -n 's/.*featureCounts v//p' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.txt
    touch ${prefix}.featureCounts.txt.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: 2.1.1
    END_VERSIONS
    """
}
