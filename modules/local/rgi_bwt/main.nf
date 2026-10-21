process RGI_BWT {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.8--pyh05cac1d_0' :
        'quay.io/biocontainers/rgi:6.0.8--pyh05cac1d_0' }"
    
    label 'process_medium'
    
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)
        path(card_db)

    output:
        tuple val(meta), path("*.allele_mapping_data.txt"), emit: allele_mapping
        tuple val(meta), path("*.gene_mapping_data.txt"), emit: gene_mapping
        tuple val(meta), path("*.sorted.length_100.bam"), emit: bam, optional: true
        path("*.overall_mapping_stats.txt"), emit: overall_stats
        path("*.reference_mapping_stats.txt"), emit: reference_stats
        path("*.artifacts_mapping_stats.txt"), emit: artifacts_stats, optional: true
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"
        def aligner = params.rgi_aligner ?: 'kma'
        def wildcard = params.rgi_include_wildcard.toString().toBoolean() ? '--include_wildcard' : ''
        
        // rgi bwt only accepts --read_one/--read_two, so read-level RGI runs on
        // the paired reads only; singleton reads (reads[2], carried since the
        // audit #9 QC rework) are deliberately excluded (documented in the manual)
        def read_one = reads[0]
        def read_two = reads.size() > 1 ? reads[1] : ''
        
        """
        set +u
        if [ -f /usr/local/env-execute ]; then
            # env-execute resets PATH; keep the pipeline's bin/ on it
            pipeline_path="\$PATH"
            source /usr/local/env-execute
            export PATH="\$PATH:\$pipeline_path"
        fi
        set -u

        # RGI --local looks for localDB/ in the current directory. Symlink the
        # DB contents per file instead of copying the multi-GB CARD/WildCARD DB
        # per task (audit #17): cp -rs makes real directories and a symlink per
        # file, so reads go through the symlinks, while any file RGI creates at
        # runtime (an aligner index under localDB/bwt/) lands in this task's
        # localDB/, not the shared published database.
        if [ ! -d ${card_db}/localDB ]; then
            echo "ERROR: ${card_db} does not contain a localDB directory (expected an 'rgi load --local' database)" >&2
            exit 1
        fi
        cp -rs "\$(readlink -f ${card_db})/localDB" localDB
        
        # Run RGI bwt through bin/rgi_bwt_kma_fix.py: RGI writes every read
        # the aligner reports, unplaced ones included (over 99% of a
        # metagenome), to the SAM and loads the read table into memory (~390
        # bytes per clean read on disk, ~1 KB in memory); the wrapper keeps
        # the reads placed on a reference (for KMA grouped by name, so samtools
        # fixmate pairs every mate pair). Remove it when RGI filters the SAM
        # itself. --clean deletes the *.temp.* files
        # when rgi finishes; the files RGI_KMER and RGI_REPORT read have no
        # .temp in their names and are kept.
        rgi_bwt_kma_fix.py bwt \\
            --read_one ${read_one} \\
            ${read_two ? "--read_two ${read_two}" : ""} \\
            --aligner ${aligner} \\
            --threads ${task.cpus} \\
            --output_file ${prefix}_rgi_bwt \\
            --local \\
            --clean \\
            ${wildcard} \\
            ${args} \\
            2>${prefix}_rgi_bwt.log

        # Low-coverage fallback (audit #8): if RGI produced no mapping results,
        # create every declared non-optional output (NA allele/gene tables plus
        # both stats files). No BAM is created - the bam emit is optional and
        # RGI_KMER is skipped for this sample (no alignments to query).
        if [ ! -f "${prefix}_rgi_bwt.allele_mapping_data.txt" ]; then
            echo "WARNING: RGI bwt did not produce expected output files - creating dummy outputs for low coverage sample"
            rgi_bwt_lowcov_fallback.sh ${prefix}_rgi_bwt
        fi

        # Create versions file (rgi --version prints the usage; KMA has no
        # --version and bwa none at all: bare bwa exits 1, hence || true
        # under pipefail)
        case ${aligner} in
            kma)  aligner_version=\$(kma -v 2>&1 | sed -n 's/^KMA-//p' || true) ;;
            bwa)  aligner_version=\$(bwa 2>&1 | sed -n 's/^Version: //p' || true) ;;
            *)    aligner_version=\$(bowtie2 --version 2>/dev/null | sed -n '1s/.* version //p' || true) ;;
        esac
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: \$(rgi main --version 2>/dev/null || echo "unknown")
            ${aligner}: \${aligner_version:-unknown}
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"
        """
        touch ${prefix}_rgi_bwt.allele_mapping_data.txt
        touch ${prefix}_rgi_bwt.gene_mapping_data.txt
        touch ${prefix}_rgi_bwt.sorted.length_100.bam
        touch ${prefix}_rgi_bwt.overall_mapping_stats.txt
        touch ${prefix}_rgi_bwt.reference_mapping_stats.txt
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.8
            kma: 1.6.11
        END_VERSIONS
        """
}
