process RGI_LOAD_WILDCARD {
    container 'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0'
    
    label 'process_medium'

    input:
        path(card_db)
        path(wildcard_dir)

    output:
        path("card_db_with_wildcard"), emit: card_db
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        # Copy CARD database to working directory
        cp -r ${card_db} card_db_with_wildcard
        cd card_db_with_wildcard

        # Copy WildCARD directory
        cp -r ${wildcard_dir} wildcard

        # Get CARD version from JSON
        CARD_VERSION=\$(grep '"version"' card.json | head -1 | sed 's/.*: "\\(.*\\)".*/\\1/')
        
        # Get CARD FASTA filename
        CARD_FASTA=\$(ls card_database_v*.fasta | head -1)
        
        if [ -z "\$CARD_FASTA" ]; then
            echo "ERROR: CARD annotation file not found in provided database"
            exit 1
        fi

        # Run wildcard annotation
        echo "Processing WildCARD annotations (this may take 30-60 minutes)..."
        rgi wildcard_annotation \\
            -i wildcard \\
            --card_json card.json \\
            -v \$CARD_VERSION \\
            > wildcard_annotation.log 2>&1
        
        # Get wildcard FASTA filename
        WILDCARD_FASTA=\$(ls wildcard_database_v*.fasta | head -1)
        
        if [ -z "\$WILDCARD_FASTA" ]; then
            echo "ERROR: WildCARD annotation file not created"
            exit 1
        fi
        
        # Load WildCARD data into existing CARD database
        rgi load \\
            --card_json card.json \\
            --wildcard_annotation \$WILDCARD_FASTA \\
            --wildcard_index wildcard/index-for-model-sequences.txt \\
            --card_annotation \$CARD_FASTA \\
            --local

        # Verify database is loaded
        rgi database --version --local > database_version.txt

        # Create versions file
        cat <<-END_VERSIONS > ../versions.yml
        "${task.process}":
            rgi: \$(rgi --version 2>&1 | grep -oP 'rgi \\K[0-9.]+' || echo "unknown")
            card: \$(grep '"version"' card.json | head -1 | sed 's/.*: "\\(.*\\)".*/\\1/')
        END_VERSIONS

        cd ..
        """

    stub:
        """
        mkdir -p card_db_with_wildcard
        touch card_db_with_wildcard/card.json
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.3
            card: 3.2.9
        END_VERSIONS
        """
}
