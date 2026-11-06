/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIN_DEPTH Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAG depth per sample from CALCULATE_DEPTH's jgi_summarize_bam_contig_depths
    table, the depth MetaBAT2 bins on (alignments at >= 97 % identity, contig
    ends trimmed): per refined bin, the contig-length-weighted mean of its
    contigs' mean depths (roadmap Phase 15, B1). One table per sample column:
    a per-sample assembly has one, a co-assembly one per sample mapped to it.

    Input:
        tuple val(meta), path(depth), path(bins), val(samples)
            depth:   <id>_depth.txt (contigName, contigLen, totalAvgDepth, then
                     <sample>_all_reads and its -var column per BAM)
            bins:    refined bin directory (*.fa)
            samples: the sample ids the depth columns must match, in any order

    Output:
        bin_depth: tuple val(meta), path("*_bin_depth.tsv") - sample, bin_id,
                   Total_contigs, average_cov
        versions:  path(versions.yml)
----------------------------------------------------------------------------------------
*/

process BIN_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabat2:2.18--h38e344b_2' :
        'quay.io/biocontainers/metabat2:2.18--h38e344b_2' }"

    input:
    tuple val(meta), path(depth), path(bins), val(samples)

    output:
    tuple val(meta), path("*_bin_depth.tsv"), emit: bin_depth
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # Bins in glob order, bin id = FASTA basename, contig = first word of the
    # header. A contig already placed in an earlier bin is skipped with a
    # warning (refined bins should never share contigs), so each contig counts
    # once; a bin left without contigs reports 0 contigs / 0 depth
    bin_files=()
    for bin in ${bins}/*.fa; do
        [ -e "\${bin}" ] && bin_files+=("\${bin}")
    done
    awk -F '\\t' -v expected="${samples.join(' ')}" '
        BEGIN {
            for (i = 1; i < ARGC - 1; i++) {
                file = ARGV[i]
                bin = file
                sub(/.*\\//, "", bin)
                sub(/\\.fa\$/, "", bin)
                order[++n_bins] = bin
                n_contigs[bin] = 0
                while ((getline line < file) > 0) {
                    if (substr(line, 1, 1) != ">") continue
                    split(substr(line, 2), words, /[ \\t]/)
                    name = words[1]
                    if (name in contig_bin) {
                        print "WARNING: contig " name " of bin " bin " is already in bin " contig_bin[name] " - skipped" > "/dev/stderr"
                        continue
                    }
                    contig_bin[name] = bin
                    n_wanted++
                }
                close(file)
                ARGV[i] = ""
            }
        }
        FNR == 1 {
            # depth columns: 4, 6, ... named after the BAM (<sample>_all_reads)
            n_samples = 0
            for (c = 4; c <= NF; c += 2) {
                s = \$c
                sub(/\\.bam\$/, "", s)
                sub(/_all_reads\$/, "", s)
                sample[++n_samples] = s
                col[n_samples] = c
            }
            n_expected = split(expected, want, " ")
            for (k = 1; k <= n_expected; k++) wanted[want[k]] = 1
            for (k = 1; k <= n_samples; k++) {
                if (!(sample[k] in wanted)) { print "ERROR: depth column " \$col[k] " is not one of the samples " expected > "/dev/stderr"; failed = 1; exit 1 }
                seen_sample[sample[k]] = 1
            }
            for (k = 1; k <= n_expected; k++) {
                if (!(want[k] in seen_sample)) { print "ERROR: no depth column for sample " want[k] > "/dev/stderr"; failed = 1; exit 1 }
            }
            next
        }
        \$1 in contig_bin {
            bin = contig_bin[\$1]
            n_contigs[bin]++
            found++
            len_sum[bin] += \$2
            for (k = 1; k <= n_samples; k++) depth_sum[bin, k] += \$2 * \$(col[k])
        }
        END {
            if (failed) exit 1
            if (found != n_wanted) { print "ERROR: " n_wanted - found " binned contig(s) are not in the depth table" > "/dev/stderr"; exit 1 }
            for (k = 1; k <= n_samples; k++) {
                out = sample[k] "_bin_depth.tsv"
                printf "sample\\tbin_id\\tTotal_contigs\\taverage_cov\\n" > out
                for (i = 1; i <= n_bins; i++) {
                    bin = order[i]
                    mean = len_sum[bin] > 0 ? depth_sum[bin, k] / len_sum[bin] : 0
                    # %s formats numbers as print does (%.6g, integers as such)
                    printf "%s\\t%s\\t%s\\t%s\\n", sample[k], bin, n_contigs[bin], mean > out
                }
                close(out)
            }
        }
    ' "\${bin_files[@]}" ${depth}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: \$(busybox 2>&1 | head -1 | sed -E 's/BusyBox v([0-9.]+).*/\\1/')
    END_VERSIONS
    """

    stub:
    """
    for s in ${samples.join(' ')}; do
        printf "sample\\tbin_id\\tTotal_contigs\\taverage_cov\\n" > "\${s}_bin_depth.tsv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: 1.36.1
    END_VERSIONS
    """
}
