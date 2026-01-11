process QFILTER {

    label 'process_single'

    container 'quay.io/ffuentessantander/r_reports:1.1'

    input:
        tuple val(meta), path(reads), path(json)

    // Output channels for this process

    output:
        tuple val(meta), path(reads), path("after_reads_fr.txt"), emit: qfilter
	path("*_fastp_report.tsv"), emit: reads_report

    script:
        def prefix = "${meta.id}"

        """
        # Handle both nf-core naming (*.fastp.json) and legacy naming (*_report.json)
        json_file=\$(ls *.json 2>/dev/null | head -1)
        
        cat \$json_file | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "after_reads.+" | cut -d':' -f2 | tr -d '\\n' > after_reads_fr.txt
        cat \$json_file | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "before.+" | cut -d':' -f2 | tr -d '\\n' > before_reads_fr.txt
        after_reads_fr=\$(cat after_reads_fr.txt)
        before_reads_fr=\$(cat before_reads_fr.txt)

        # Check for singleton JSON (legacy format)
        singleton_json=\$(ls *Singleton*.json 2>/dev/null | head -1 || echo "")
        
        if [[ -n "\$singleton_json" && -s "\$singleton_json" ]]; then
            cat \$singleton_json | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "after_reads.+" | cut -d':' -f2 | tr -d '\\n' > after_reads_Singleton
            cat \$singleton_json | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "before_reads.+" | cut -d':' -f2 | tr -d '\\n' > before_reads_Singleton
            after_reads_singleton=\$(cat after_reads_Singleton)
            before_reads_singleton=\$(cat before_reads_Singleton)
            
            echo -e "Id\\tRaw reads\\tFastp\\tRaw singletons\\tFastp singletons" > ${prefix}_fastp_report.tsv
            echo -e "${prefix}\\t\${before_reads_fr}\\t\${after_reads_fr}\\t\${before_reads_singleton}\\t\${after_reads_singleton}" >> ${prefix}_fastp_report.tsv 
        else
            echo -e "Id\\tRaw reads\\tFastp" > ${prefix}_fastp_report.tsv
            echo -e "${prefix}\\t\${before_reads_fr}\\t\${after_reads_fr}" >> ${prefix}_fastp_report.tsv
        fi
        """ 
}
