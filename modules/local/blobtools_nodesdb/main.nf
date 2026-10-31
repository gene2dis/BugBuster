process BLOBTOOLS_NODESDB {
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/blobtools:1.1.1--py_1' :
        'quay.io/biocontainers/blobtools:1.1.1--py_1' }"

    tag "nodesdb"
    label 'process_single'

    input:
        path(tax_files)

    output:
        path("nodesDB.txt"), emit: nodesdb
        path "versions.yml", emit: versions

    script:
        """
        # tax_files: nodes.dmp and names.dmp of the Metabuli database, whose
        # taxids the hits files carry
        nodes_dmp=\$(find -L . -maxdepth 2 -name nodes.dmp | head -n 1)
        names_dmp=\$(find -L . -maxdepth 2 -name names.dmp | head -n 1)

        # 2025 NCBI taxdumps call the top rank 'domain' (and Viruses an
        # 'acellular root'): renamed 'superkingdom' in a task-dir copy, or
        # BlobTools writes 'undef' as every superkingdom
        blobtools_nodes_fix.sh \${nodes_dmp} > nodes_blobtools.dmp

        # The parse 'blobtools create --nodes --names --db' does, once for
        # every sample: 'blobtools nodesdb' would write into the blobtools
        # package, read-only in a singularity/apptainer image (design doc Q21)
        python -c "import blobtools.BtIO as BtIO; BtIO.parseNodesDB(nodes='nodes_blobtools.dmp', names='\${names_dmp}', nodesDB='nodesDB.txt', nodesDBdefault=None)"
        rm nodes_blobtools.dmp

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: \$(blobtools --version 2>&1 | tail -1 || echo 1.1.1)
        END_VERSIONS
        """

    stub:
        """
        printf '# nodes_count = 1\\n1\\tno rank\\troot\\t1\\n' > nodesDB.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            blobtools: 1.1.1
        END_VERSIONS
        """
}
