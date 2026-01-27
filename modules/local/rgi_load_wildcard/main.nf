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
        set +u
        if [ -f /usr/local/env-execute ]; then
            source /usr/local/env-execute
        fi
        set -u

        # Copy CARD database to working directory (avoid circular symlinks)
        mkdir -p card_db_with_wildcard
        
        # Use rsync to safely copy and dereference symlinks while avoiding loops
        if command -v rsync &> /dev/null; then
            # Exclude the circular 'rgi' symlink that points to parent directory
            rsync -rL --copy-unsafe-links --exclude='rgi' ${card_db}/ card_db_with_wildcard/
        else
            # Fallback: copy without dereferencing symlinks to avoid loops
            cp -r ${card_db}/* card_db_with_wildcard/
            # Remove any problematic self-referencing symlinks
            find card_db_with_wildcard -type l -exec sh -c 'readlink "\$1" | grep -q "^\$(dirname "\$1")" && rm "\$1"' _ {} \\;
        fi
        
        # Additional cleanup: remove any circular symlinks that may have been copied
        find card_db_with_wildcard -type l -name 'rgi' -delete 2>/dev/null || true
        
        # Copy WildCARD directory contents before changing directory
        mkdir -p card_db_with_wildcard/wildcard_data
        if [ -d "${wildcard_dir}" ]; then
            if command -v rsync &> /dev/null; then
                rsync -rL --copy-unsafe-links ${wildcard_dir}/ card_db_with_wildcard/wildcard_data/
            else
                cp -r ${wildcard_dir}/* card_db_with_wildcard/wildcard_data/ 2>/dev/null || true
            fi
        else
            echo "ERROR: WildCARD directory not found: ${wildcard_dir}"
            exit 1
        fi
        
        cd card_db_with_wildcard

        # Get CARD version from JSON using Python (more reliable than grep/sed)
        CARD_VERSION=\$(python3 -c "import json; data=json.load(open('card.json')); print(data.get('_version', ''))" 2>/dev/null || echo "")
        
        # Get CARD FASTA filename
        CARD_FASTA=\$(ls card_database_v*.fasta | head -1)
        
        if [ -z "\$CARD_FASTA" ]; then
            echo "ERROR: CARD annotation file not found in provided database"
            exit 1
        fi

        # Check if wildcard database already exists
        WILDCARD_FASTA=\$(ls wildcard_database_v*.fasta 2>/dev/null | head -1)
        
        if [ -z "\$WILDCARD_FASTA" ]; then
            # Run wildcard annotation only if not already present
            echo "Processing WildCARD annotations (this may take 30-60 minutes)..."
            rgi wildcard_annotation \\
                -i wildcard_data \\
                --card_json card.json \\
                -v \$CARD_VERSION \\
                > wildcard_annotation.log 2>&1
            
            # Get wildcard FASTA filename after creation
            WILDCARD_FASTA=\$(ls wildcard_database_v*.fasta | head -1)
            
            if [ -z "\$WILDCARD_FASTA" ]; then
                echo "ERROR: WildCARD annotation file not created"
                cat wildcard_annotation.log
                exit 1
            fi
        else
            echo "WildCARD database already exists: \$WILDCARD_FASTA"
        fi
        
        # Load WildCARD data into existing CARD database
        rgi load \\
            --card_json card.json \\
            --wildcard_annotation \$WILDCARD_FASTA \\
            --wildcard_index wildcard_data/index-for-model-sequences.txt \\
            --card_annotation \$CARD_FASTA \\
            --local

        # Verify database is loaded
        rgi database --version --local > database_version.txt

        # Create versions file
        cat <<-END_VERSIONS > ../versions.yml
        "${task.process}":
            rgi: \$(rgi --version 2>&1 | sed -n 's/.*rgi \\([0-9.]*\\).*/\\1/p' || echo "unknown")
            card: \$(python3 -c "import json; data=json.load(open('card.json')); print(data.get('_version', 'unknown'))" 2>/dev/null || echo "unknown")
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
