process BEDTOOLS {
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    tag "${meta.id}"
    label 'process_single'

    input:
        tuple val(meta), path(bam_files)

    // Output channels for this process

    output:
        path("*_bin_depth.tsv"), emit: bin_depth
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        printf "sample\\tbin_id\\tTotal_contigs\\taverage_cov\\n" > ${prefix}_bin_depth.tsv
        for bam in *.bam; do
            # BAMs are named <sample>_<bin>_all_reads.bam by
            # BOWTIE2_SAMTOOLS_DEPTH; strip prefix and suffix instead of the
            # old greedy-sed parse that only worked for bins named bin.N
            bin_name=\${bam#${prefix}_}
            bin_name=\${bin_name%_all_reads.bam}
            genomeCoverageBed -ibam \${bam} > \${bin_name}_cov.tsv
            # Skip the whole-genome histogram rows (chrom = "genome"): they are
            # a summary of the contig rows, not another contig
            cat \${bin_name}_cov.tsv | awk '\$1 != "genome" {print \$1, \$2*\$3, \$4, \$5}' | awk '{sum[\$1] += \$2; contig_length[\$1] = \$3} END { for (nombre in sum) { promedio = sum[nombre] / contig_length[nombre]; print nombre, promedio, sum[nombre], contig_length[nombre]; } }' | awk '{sum += \$2; count++; } END { if (count > 0) { promedio = sum / count } else { count = 0; promedio = 0 }; print "tmp1\ttmp2\t",count,"\t",promedio; }' > bin_depth_tmp.tsv
            cat bin_depth_tmp.tsv | sed "s/tmp1/${prefix}/g" | sed "s/tmp2/\${bin_name}/g" >> ${prefix}_bin_depth.tsv
        done

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

/*
    BEDTOOLS_BIN_DEPTH (--read_aligner bwamem3): same per-bin coverage table
    as BEDTOOLS, from ONE BAM of the sample's reads mapped against all bins
    at once (BINS_REFERENCE + BWAMEM3_SAMTOOLS). contig_to_bin.tsv assigns
    each contig to its bin; a bin row with an empty contig field is a bin
    without sequences and reports 0 contigs / 0 coverage, like an empty
    per-bin BAM does in BEDTOOLS. Per-bin values keep BEDTOOLS' semantics:
    the mean over the bin's contigs of each contig's mean depth.
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
        printf "sample\\tbin_id\\tTotal_contigs\\taverage_cov\\n" > ${prefix}_bin_depth.tsv
        genomeCoverageBed -ibam ${bam} > ${prefix}_cov.tsv
        # Skip the whole-genome histogram rows (chrom = "genome"): they are a
        # summary of the contig rows, not another contig
        awk -F '\\t' -v sample="${prefix}" '
            FILENAME == ARGV[1] {
                if (!(\$2 in n_contigs)) { order[++n_bins] = \$2; n_contigs[\$2] = 0 }
                if (\$1 != "") contig_bin[\$1] = \$2
                next
            }
            \$1 != "genome" { depth_sum[\$1] += \$2 * \$3; contig_length[\$1] = \$4 }
            END {
                for (contig in contig_length) {
                    if (contig in contig_bin) {
                        bin = contig_bin[contig]
                        cov_sum[bin] += depth_sum[contig] / contig_length[contig]
                        n_contigs[bin]++
                    }
                }
                for (i = 1; i <= n_bins; i++) {
                    bin = order[i]
                    promedio = n_contigs[bin] > 0 ? cov_sum[bin] / n_contigs[bin] : 0
                    print sample "\\t" bin "\\t", n_contigs[bin], "\\t", promedio
                }
            }
        ' ${contig_to_bin} ${prefix}_cov.tsv >> ${prefix}_bin_depth.tsv
        rm -f ${prefix}_cov.tsv

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
