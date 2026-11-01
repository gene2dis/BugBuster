/*
    BIN_CONTIG_TABLE: contig -> bin table of the refined bins for
    BEDTOOLS_BIN_DEPTH (co-assembly bin depth).
*/
process BIN_CONTIG_TABLE {
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    tag "bins"
    label 'process_single'

    input:
        path(bins)

    output:
        path("contig_to_bin.tsv"), emit: table
        path "versions.yml", emit: versions

    script:
        """
        set -euo pipefail

        # contig_to_bin.tsv: <contig> TAB <bin>, bins in glob order, contig =
        # first word of the header, bin id = FASTA basename. A contig already
        # placed in an earlier bin is skipped with a warning (refined bins
        # should never share contigs), so each contig counts once. A bin that
        # adds no contig (no sequences, or only contigs of earlier bins) gets
        # one row with an empty contig field, so it still reports a
        # 0-coverage row downstream.
        : > contig_to_bin.tsv
        bin_files=()
        for bin in ${bins}/*.fa; do
            [ -e "\${bin}" ] && bin_files+=("\${bin}")
        done
        awk '
            BEGIN {
                for (i = 1; i < ARGC; i++) {
                    file = ARGV[i]
                    bin = file
                    sub(/.*\\//, "", bin)
                    sub(/\\.fa\$/, "", bin)
                    kept = 0
                    while ((getline line < file) > 0) {
                        if (substr(line, 1, 1) != ">") continue
                        split(substr(line, 2), words, /[ \\t]/)
                        name = words[1]
                        if (name in seen) {
                            print "WARNING: contig " name " of bin " bin " is already in another bin - skipped" > "/dev/stderr"
                            continue
                        }
                        seen[name] = 1
                        kept++
                        print name "\\t" bin > "contig_to_bin.tsv"
                    }
                    close(file)
                    if (!kept) print "\\t" bin > "contig_to_bin.tsv"
                }
                exit
            }
        ' "\${bin_files[@]}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            busybox: \$(busybox 2>&1 | head -1 | sed -E 's/BusyBox v([0-9.]+).*/\\1/')
        END_VERSIONS
        """

    stub:
        """
        printf 'contig_1\\tbin.1\\n' > contig_to_bin.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            busybox: 1.36.1
        END_VERSIONS
        """
}

/*
    BEDTOOLS_BIN_DEPTH: co-assembly bin coverage per sample, from the
    sample's alignment against the whole co-assembly (the same BAM that
    counts genes and gives MetaBAT2 its depths; issue #28). Only the contigs
    in contig_to_bin.tsv (BIN_CONTIG_TABLE) count, each in its bin; a bin row
    with an empty contig field is a bin without sequences and reports
    0 contigs / 0 coverage. Per bin: the mean over its contigs of each
    contig's mean depth.
*/
process BEDTOOLS_BIN_DEPTH {
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    tag "${meta.id}"
    label 'process_single'

    input:
        tuple val(meta), path(bam), path(contig_to_bin)

    output:
        path("*_bin_depth.tsv"), emit: bin_depth
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        set -euo pipefail

        printf "sample\\tbin_id\\tTotal_contigs\\taverage_cov\\n" > ${prefix}_bin_depth.tsv
        # The histogram streams into awk (the co-assembly BAM's header holds
        # every contig) and only binned contigs are kept, which also skips
        # the whole-genome summary rows (chrom = "genome")
        genomeCoverageBed -ibam ${bam} | awk -F '\\t' -v sample="${prefix}" '
            FILENAME == ARGV[1] {
                if (!(\$2 in n_contigs)) { order[++n_bins] = \$2; n_contigs[\$2] = 0 }
                if (\$1 != "") contig_bin[\$1] = \$2
                next
            }
            \$1 in contig_bin { depth_sum[\$1] += \$2 * \$3; contig_length[\$1] = \$4 }
            END {
                for (contig in contig_length) {
                    bin = contig_bin[contig]
                    cov_sum[bin] += depth_sum[contig] / contig_length[contig]
                    n_contigs[bin]++
                }
                for (i = 1; i <= n_bins; i++) {
                    bin = order[i]
                    promedio = n_contigs[bin] > 0 ? cov_sum[bin] / n_contigs[bin] : 0
                    # %s formats numbers as print did (%.6g, integers as
                    # such); print with commas put a space around each field
                    printf "%s\\t%s\\t%s\\t%s\\n", sample, bin, n_contigs[bin], promedio
                }
            }
        ' ${contig_to_bin} - >> ${prefix}_bin_depth.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed 's/bedtools v//')
        END_VERSIONS
	"""

    stub:
        def prefix = "${meta.id}"

        """
        printf "sample\\tbin_id\\tTotal_contigs\\taverage_cov\\n" > ${prefix}_bin_depth.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: 2.31.1
        END_VERSIONS
        """
}
