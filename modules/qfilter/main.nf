process QFILTER {

    label 'process_single'

    container 'quay.io/ffuentessantander/r_reports:1.1'

    input:
        tuple val(meta), path(reads), path(json)

    // Aqui creo los canales de salida de este proceso.

    output:
        tuple val(meta), path(reads), path("after_reads_fr.txt"), emit: qfilter
	path("*_fastp_report.tsv"), emit: reads_report

    script:
        def prefix = "${meta.id}"

        """

        cat ${prefix}_report.json | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "after_reads.+" | cut -d':' -f2 | tr -d '\\n' > after_reads_fr.txt
	cat ${prefix}_report.json | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "before.+" | cut -d':' -f2 | tr -d '\\n' > before_reads_fr.txt
        after_reads_fr=`cat after_reads_fr.txt`
        before_reads_fr=`cat before_reads_fr.txt`

	if [[ `cat ${prefix}_Singleton_report.json` != "" ]]; then
		cat ${prefix}_Singleton_report.json | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "after_reads.+" | cut -d':' -f2 | tr -d '\\n' > after_reads_Singleton
		cat ${prefix}_Singleton_report.json | grep -zoP '.+_filtering.+\\n.+'| sed 's/"//g' | awk '{print \$1}' | grep -Pa 'total_reads:[0-9].+' | sed -e 's/,//g' -e '1s/total/before/' -e '2s/total/after/g'| tr '\\0' '\\n' | grep -Poa "before_reads.+" | cut -d':' -f2 | tr -d '\\n' > before_reads_Singleton
		after_reads_singleton=`cat after_reads_Singleton`
		before_reads_singleton=`cat before_reads_Singleton`
	fi

	if [[ `cat ${prefix}_Singleton_report.json` != "" ]]; then
		echo "Id\tRaw reads\tFastp\tRaw singletons\tFastp singletons" > ${prefix}_fastp_report.tsv
		echo "${prefix}\t\${before_reads_fr}\t\${after_reads_fr}\t\${before_reads_singleton}\t\${after_reads_singleton}" >> ${prefix}_fastp_report.tsv 

	else
		echo "Id\tRaw reads\tFastp" > ${prefix}_fastp_report.tsv
		echo "${prefix}\t\${before_reads_fr}\t\${after_reads_fr}" >> ${prefix}_fastp_report.tsv
	fi
	
	""" 
}
