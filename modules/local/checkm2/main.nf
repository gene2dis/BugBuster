/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKM2 Quality Assessment Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Bin quality assessment using CheckM2
    
    Evaluates completeness and contamination of bins from multiple binners
    
    Input:
        tuple val(meta), path(metabat2), path(semibin), path(comebin), path(metawrap), path(checkm_db)
    
    Output:
        all_reports: tuple val(meta), path(quality_reports)
        metawrap_report: tuple val(meta), path(metawrap_quality_report)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CHECKM2 {
    tag "${meta.id}"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.1.0--pyh7e72e81_1' :
        'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1' }"

    input:
    tuple val(meta), path(metabat2), path(semibin), path(comebin), path(metawrap), path(checkm_db)

    output:
    tuple val(meta), path("*quality_report.tsv"), emit: all_reports
    tuple val(meta), path("*_metawrap_quality_report.tsv"), emit: metawrap_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = meta.id
    
    // Input validation
    if (!meta || !meta.id) {
        error "CHECKM2: meta.id is required"
    }
    if (!checkm_db) {
        error "CHECKM2: CheckM2 database path is required"
    }
    
    """
    set -euo pipefail
    
    # Run CheckM2 on each binner's output
    for binner in metabat semibin comebin metawrap; do
        case \$binner in
            metabat)
                bin_dir="${metabat2}"
                ;;
            semibin)
                bin_dir="${semibin}"
                ;;
            comebin)
                bin_dir="${comebin}"
                ;;
            metawrap)
                bin_dir="${metawrap}"
                ;;
        esac
        
        checkm2 predict \\
            --threads ${task.cpus} \\
            --input \$bin_dir \\
            --output-directory ${prefix}_\${binner}_checkm2_report \\
            --database_path ${checkm_db} \\
            -x .fa
        
        mv ${prefix}_\${binner}_checkm2_report/quality_report.tsv ${prefix}_\${binner}_quality_report.tsv
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version 2>&1 | sed 's/checkm2: version //')
    END_VERSIONS
    """
    
    stub:
    prefix = meta.id
    """
    # Create stub quality report files
    for binner in metabat semibin comebin metawrap; do
        cat > ${prefix}_\${binner}_quality_report.tsv << 'EOF'
Name	Completeness	Contamination
bin.1	95.5	2.1
bin.2	87.3	1.5
EOF
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: 1.1.0
    END_VERSIONS
    """
}
