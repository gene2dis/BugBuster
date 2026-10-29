process BLOBTOOLS {
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/blobtools:1.1.1--py_1' :
        'quay.io/biocontainers/blobtools:1.1.1--py_1' }"

    tag "${meta.id}"
    label 'process_single'

    input:
        tuple val(meta), path(contigs), path(hits), path(bam), path(bam_bai)
        path(tax_files)

    output:
        tuple val(meta), path("*_Blob_tabl*"), emit: blob_table
        path("*_Blob_tabl*"), emit: only_blob
        path "versions.yml", emit: versions

    script:
        def prefix = "${meta.id}"

        """
        # tax_files: nodes.dmp and names.dmp of the Metabuli database, whose
        # taxids the hits file carries
        nodes_dmp=\$(find -L . -maxdepth 2 -name nodes.dmp | head -n 1)
        names_dmp=\$(find -L . -maxdepth 2 -name names.dmp | head -n 1)

        # 2025 NCBI taxdumps call the top rank 'domain' (and Viruses an
        # 'acellular root'): renamed 'superkingdom' in a task-dir copy, or
        # BlobTools writes 'undef' as every superkingdom
        blobtools_nodes_fix.sh \${nodes_dmp} > nodes_blobtools.dmp

        # --db: write the parsed nodesDB into the task dir; the default is
        # inside the blobtools package, which is read-only in a
        # singularity/apptainer image (design doc Q21)
        blobtools create \\
                  --infile ${contigs} \\
                  --hitsfile ${hits} \\
                  --nodes nodes_blobtools.dmp \\
                  --names \${names_dmp} \\
                  --db nodesDB.txt \\
                  --bam ${bam} \\
                  --out ${prefix}

        blobtools view \\
                  --input ${prefix}.blobDB.json \\
                  --out ${prefix}_Blob_table \\
                  --rank all

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: \$(blobtools --version 2>&1 | tail -1 || echo 1.1.1)
        END_VERSIONS
        """

    stub:
        def prefix = "${meta.id}"

        """
        touch ${prefix}_Blob_table.blobDB.table.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: 1.1.1
        END_VERSIONS
        """
}
