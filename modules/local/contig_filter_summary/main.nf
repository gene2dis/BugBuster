/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONTIG_FILTER_SUMMARY Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Aggregate contig filtering reports to track samples with empty/filtered contigs
    
    Generates a comprehensive summary report showing:
    - Samples where all contigs were filtered out by BBMAP
    - Samples where empty contigs were detected by BOWTIE2_SAMTOOLS
    - Filtering statistics (contigs before/after, thresholds)
    
    Input:
        path(filter_reports) - Collection of BBMAP filter reports
        path(empty_reports)  - Collection of BOWTIE2_SAMTOOLS empty contig reports
    
    Output:
        summary_report: path(contig_filtering_summary.txt)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CONTIG_FILTER_SUMMARY {
    label 'process_single'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"
    
    publishDir "${params.output}/pipeline_info", mode: params.publish_dir_mode

    input:
    path(filter_reports)
    path(empty_reports)

    output:
    path "contig_filtering_summary.txt", emit: summary_report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail
    
    # Create summary report header
    cat > contig_filtering_summary.txt <<-HEADER
# Contig Filtering Summary Report
# Generated: \$(date -u +"%Y-%m-%d %H:%M:%S UTC")
# Pipeline: BugBuster
#
# This report tracks samples where contig filtering resulted in empty outputs
#
================================================================================

HEADER
    
    # Section 1: BBMAP Filtering Statistics
    echo "## BBMAP Contig Filtering Statistics" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    
    if [ -f "${filter_reports}" ] || ls *_filter_report.txt 1> /dev/null 2>&1; then
        printf "Sample_ID\\tContigs_Before\\tContigs_After\\tMin_Length\\tAll_Filtered\\n" >> contig_filtering_summary.txt
        echo "------------------------------------------------------------------------" >> contig_filtering_summary.txt
        
        for report in *_filter_report.txt; do
            if [ -f "\$report" ]; then
                # Extract values from report
                sample=\$(grep "^sample_id" "\$report" | cut -f2)
                before=\$(grep "^contigs_before_filter" "\$report" | cut -f2)
                after=\$(grep "^contigs_after_filter" "\$report" | cut -f2)
                threshold=\$(grep "^min_length_threshold" "\$report" | cut -f2)
                filtered=\$(grep "^all_contigs_filtered" "\$report" | cut -f2)
                
                printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\$sample" "\$before" "\$after" "\$threshold" "\$filtered" >> contig_filtering_summary.txt
            fi
        done
    else
        echo "No BBMAP filter reports found" >> contig_filtering_summary.txt
    fi
    
    echo "" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    
    # Section 2: Samples with Empty Contigs (caught by BOWTIE2_SAMTOOLS)
    echo "## Samples with Empty Contigs Detected at Alignment Stage" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    
    if [ -f "${empty_reports}" ] || ls *_empty_contig_report.txt 1> /dev/null 2>&1; then
        printf "Sample_ID\\tStage\\tTimestamp\\tAction_Taken\\n" >> contig_filtering_summary.txt
        echo "------------------------------------------------------------------------" >> contig_filtering_summary.txt
        
        for report in *_empty_contig_report.txt; do
            if [ -f "\$report" ]; then
                sample=\$(grep "^sample_id" "\$report" | cut -f2)
                stage=\$(grep "^stage" "\$report" | cut -f2)
                timestamp=\$(grep "^timestamp" "\$report" | cut -f2-)
                action=\$(grep "^action_taken" "\$report" | cut -f2-)
                
                printf "%s\\t%s\\t%s\\t%s\\n" "\$sample" "\$stage" "\$timestamp" "\$action" >> contig_filtering_summary.txt
            fi
        done
    else
        echo "No empty contig reports found (all samples had contigs after filtering)" >> contig_filtering_summary.txt
    fi
    
    echo "" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    
    # Section 3: Summary Statistics
    echo "## Summary Statistics" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    
    total_samples=\$(ls *_filter_report.txt 2>/dev/null | wc -l || echo 0)
    empty_samples=\$(grep -l "all_contigs_filtered.*TRUE" *_filter_report.txt 2>/dev/null | wc -l || echo 0)
    
    echo "Total samples processed: \$total_samples" >> contig_filtering_summary.txt
    echo "Samples with all contigs filtered: \$empty_samples" >> contig_filtering_summary.txt
    
    if [ \$total_samples -gt 0 ]; then
        # Calculate percentage using bash arithmetic (integer only)
        percentage=\$(( (\$empty_samples * 100) / \$total_samples ))
        echo "Percentage of samples with empty contigs: \${percentage}%" >> contig_filtering_summary.txt
    fi
    
    echo "" >> contig_filtering_summary.txt
    echo "================================================================================" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    echo "Note: Samples with empty contigs may indicate:" >> contig_filtering_summary.txt
    echo "  - Low sequencing depth" >> contig_filtering_summary.txt
    echo "  - Poor read quality" >> contig_filtering_summary.txt
    echo "  - Assembly failure" >> contig_filtering_summary.txt
    echo "  - Overly strict length filtering threshold" >> contig_filtering_summary.txt
    echo "" >> contig_filtering_summary.txt
    echo "Consider reviewing these samples and adjusting --bbmap_length if needed." >> contig_filtering_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    cat > contig_filtering_summary.txt << 'EOF'
# Contig Filtering Summary Report
# Generated: 2026-01-22 21:20:00 UTC
# Pipeline: BugBuster

================================================================================

## BBMAP Contig Filtering Statistics

Sample_ID	Contigs_Before	Contigs_After	Min_Length	All_Filtered
------------------------------------------------------------------------
sample1	100	50	1000	FALSE
sample2	25	0	1000	TRUE

## Samples with Empty Contigs Detected at Alignment Stage

Sample_ID	Stage	Timestamp	Action_Taken
------------------------------------------------------------------------
sample2	BOWTIE2_SAMTOOLS	2026-01-22 21:15:00 UTC	Created empty BAM file

## Summary Statistics

Total samples processed: 2
Samples with all contigs filtered: 1
Percentage of samples with empty contigs: 50.00%

================================================================================
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: 5.1.0
    END_VERSIONS
    """
}
