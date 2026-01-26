process RGI_LOAD {
    container 'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0'
    
    label 'process_download'

    input:
        val(card_version)
        val(include_wildcard)

    output:
        path("card_db"), emit: card_db
        path("versions.yml"), emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        mkdir -p card_db
        cd card_db

        # Download CARD database
        echo "Downloading CARD database version: ${card_version}"
        if [ "${card_version}" = "latest" ]; then
            wget -O data.tar.bz2 https://card.mcmaster.ca/latest/data
        else
            wget -O data.tar.bz2 https://card.mcmaster.ca/download/0/broadstreet-v${card_version}.tar.bz2
        fi
        
        tar -xjf data.tar.bz2
        rm data.tar.bz2

        # Clean any previous RGI data
        rgi clean --local

        # Load CARD JSON
        rgi load --card_json card.json --local

        # Create annotated reference sequences
        rgi card_annotation -i card.json > card_annotation.log 2>&1

        # Get the actual filename created
        CARD_FASTA=\$(ls card_database_v*.fasta | head -1)
        
        if [ -z "\$CARD_FASTA" ]; then
            echo "ERROR: CARD annotation file not created"
            exit 1
        fi

        # Load annotated sequences
        rgi load -i card.json \\
            --card_annotation \$CARD_FASTA \\
            --local

        # Download and prepare WildCARD if requested
        if [ "${include_wildcard}" = "true" ]; then
            echo "Downloading WildCARD variants..."
            wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
            
            mkdir -p wildcard
            tar -xjf wildcard_data.tar.bz2 -C wildcard
            rm wildcard_data.tar.bz2
            
            # Decompress all files
            gunzip wildcard/*.gz || true
            
            # Get CARD version from JSON
            CARD_VERSION=\$(grep '"version"' card.json | head -1 | sed 's/.*: "\\(.*\\)".*/\\1/')
            
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
            
            # Load WildCARD data
            rgi load \\
                --card_json card.json \\
                --wildcard_annotation \$WILDCARD_FASTA \\
                --wildcard_index wildcard/index-for-model-sequences.txt \\
                --card_annotation \$CARD_FASTA \\
                --local
        fi

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
        mkdir -p card_db
        touch card_db/card.json
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rgi: 6.0.3
            card: 3.2.9
        END_VERSIONS
        """
}
