process SOURMASH {
    container 'quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0'

    maxForks 24
    label 'process_single'

    publishDir "${params.output}/workflow/${meta.id}/sourmash_${db_name}", pattern: '*.with-lineages.csv'

    input:
        tuple val(meta), path(reads), path(db), path(tax_file)
	val db_name
        val sourmash_rank

    output:
        path("*.with-lineages.csv"), emit: sourmash_gather
	path("*_report.tsv"), emit: report

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = "${meta.id}"

        """
        if [[ ${reads[2]} != null ]]; then

           case \$(echo ${sourmash_rank} | awk '{print tolower(\$0)}') in 
              genus)
                sourmash sketch dna <(zcat ${reads[0]} ${reads[1]} ${reads[2]}) -p k=21,scaled=1000,abund -o ${prefix}.sig
                sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k21_renamed.sig
                sourmash gather ${prefix}_k21_renamed.sig ${db} -k 21 -o ${prefix}_smgather_${db_name}.csv
                ;;

              species)
                sourmash sketch dna <(zcat ${reads[0]} ${reads[1]} ${reads[2]}) -p k=31,scaled=1000,abund -o ${prefix}.sig
                sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k31_renamed.sig
                sourmash gather ${prefix}_k31_renamed.sig ${db} -k 31 -o ${prefix}_smgather_${db_name}.csv
                ;;

              strain)
                sourmash sketch dna <(zcat ${reads[0]} ${reads[1]} ${reads[2]}) -p k=51,scaled=1000,abund -o ${prefix}.sig
                sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k51_renamed.sig
                sourmash gather ${prefix}_k51_renamed.sig ${db} -k 51 -o ${prefix}_smgather_${db_name}.csv
                ;;
           esac
        else

           case \$(echo ${sourmash_rank} | awk '{print tolower(\$0)}') in
              genus)
                sourmash sketch dna <(zcat ${reads[0]} ${reads[1]} ${reads[2]}) -p k=21,scaled=1000,abund -o ${prefix}.sig
                sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k21_renamed.sig
                sourmash gather ${prefix}_k21_renamed.sig ${db} -k 21 -o ${prefix}_smgather_${db_name}.csv
                ;;
 
              species)
                sourmash sketch dna <(zcat ${reads[0]} ${reads[1]}) -p k=31,scaled=1000,abund -o ${prefix}.sig
                sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k31_renamed.sig
                sourmash gather ${prefix}_k31_renamed.sig ${db} -k 31 -o ${prefix}_smgather_${db_name}.csv
              ;;

              strain)
                sourmash sketch dna <(zcat ${reads[0]} ${reads[1]}) -p k=51,scaled=1000,abund -o ${prefix}.sig
                sourmash signature rename ${prefix}.sig "${prefix}" -o ${prefix}_k51_renamed.sig
                sourmash gather ${prefix}_k51_renamed.sig ${db} -k 51 -o ${prefix}_smgather_${db_name}.csv
              ;;
           esac
        fi

        sourmash tax prepare -t ${tax_file} -o tax_file.taxonomy.sqldb -F sql       
        sourmash tax annotate -g ${prefix}_smgather_${db_name}.csv -t tax_file.taxonomy.sqldb

        classified=`awk -F, '{sum += \$5} END {print sum}' ${prefix}_smgather_${db_name}.with-lineages.csv`
        unclassified=`echo "1 \${classified}" | awk '{print \$1 - \$2}'`
	echo "Id\tSourmash DB\tUnclassified\tClassified" > ${meta.id}_${sourmash_rank}_${db_name}_report.tsv
        echo "${prefix}\t${db_name}\t\${unclassified}\t\${classified}" >> ${meta.id}_${sourmash_rank}_${db_name}_report.tsv

        rm -f *.sig
        rm -f tax_file.taxonomy.sqldb
        """
}
