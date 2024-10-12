process BEDTOOLS {
    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_1'

    publishDir "${params.output}/workflow/${meta.id}/Bins_depth", pattern: '*_bin_depth.tsv'

    label 'process_single'

    input:
        tuple val(meta), path(bam_files)

    // Aqui creo los canales de salida de este proceso.

    output:
        path("*_bin_depth.tsv"), emit: bin_depth

    script:
        def prefix = "${meta.id}"

        """ 
        echo "sample\tbin_id\tTotal_contigs\taverage_cov" > ${prefix}_bin_depth.tsv
        for bam in *.bam; do
            bin_name=`echo \${bam} | sed -E 's/.+?_b//g' | sed 's/_all_reads.bam//g' | sed 's/in/bin/g'`
            genomeCoverageBed -ibam \${bam} > \${bin_name}_cov.tsv       
            cat \${bin_name}_cov.tsv | awk '{print \$1, \$2*\$3, \$4, \$5}' | awk '{sum[\$1] += \$2; contig_length[\$1] = \$3} END { for (nombre in sum) { promedio = sum[nombre] / contig_length[nombre]; print nombre, promedio, sum[nombre], contig_length[nombre]; } }' | awk '{sum += \$2; count++; } END {promedio = sum / count; print "tmp1\ttmp2\t",count,"\t",promedio; }' > bin_depth_tmp.tsv
            cat bin_depth_tmp.tsv | sed "s/tmp1/${prefix}/g" | sed "s/tmp2/\${bin_name}/g" >> ${prefix}_bin_depth.tsv
        done 
	""" 
}
