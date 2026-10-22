/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECKM2_BATCH Quality Assessment Module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Batched bin quality assessment using CheckM2
    
    One checkm2 predict over every bin of every sample and binner (bins renamed
    <sample>__<binner>__<bin>), so the DIAMOND database is searched and the
    models are loaded once per run instead of once per binner (performance
    review E5); the report is split back into one file per sample and binner.

    Input:
        tuple val(meta_list), path(all_bins, stageAs: 'bin_?/*'), path(checkm_db)
    
    Output:
        all_reports: tuple val(meta), path(quality_reports) - per sample
        metawrap_report: tuple val(meta), path(metawrap_quality_report) - per sample (optional; only when MetaWRAP ran)
        versions: path(versions.yml)
----------------------------------------------------------------------------------------
*/

process CHECKM2_BATCH {
    tag "batch_all_samples"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.1.0--pyh7e72e81_1' :
        'quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_1' }"

    input:
    tuple val(meta_list), path(all_bins, stageAs: 'bin_?/*'), path(checkm_db)

    output:
    tuple val(meta_list), path("*_*_quality_report.tsv"), emit: all_reports
    tuple val(meta_list), path("*_metawrap_quality_report.tsv"), emit: metawrap_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.unique().join(' ')
    def metawrap_suffix = task.ext.metawrap_dir_suffix ?: "metawrap_${params.metawrap_completeness}_${params.metawrap_contamination}_bins"
    """
    set -euo pipefail

    # Every bin of every sample and binner goes into one directory as
    # <sample>__<binner>__<bin>.fa, so one checkm2 run (one DIAMOND search,
    # one model load) covers them all. bin_map.tsv maps each CheckM2 name
    # (the file name without .fa) back to its sample, binner and bin name.
    mkdir -p all_bins
    printf 'name\\tsample_id\\tbinner\\tbin\\n' > bin_map.tsv
    : > binners.txt

    # Process ALL staged bin directories - parse sample_id and binner from directory names
    # Bins are staged as bin_1/*, bin_2/*, etc. by Nextflow. Every staged dir is a
    # binner output, so identify_bin_dir.sh fails loudly on naming drift (audit #6);
    # the MetaWRAP dir name depends on params.metawrap_* and is passed in via ext.
    for staged_dir in bin_*/; do
        staged_dir=\${staged_dir%/}  # Remove trailing slash

        # Get the actual bin directory inside
        bin_dir=\$(find "\$staged_dir" -mindepth 1 -maxdepth 1 -type d -o -type l | head -1)
        [ -z "\$bin_dir" ] && continue

        id_line=\$(identify_bin_dir.sh "\$bin_dir" "${metawrap_suffix}")
        sample_id=\${id_line%%\$'\\t'*}
        binner=\${id_line##*\$'\\t'}

        echo "Processing \$bin_dir -> sample_id=\$sample_id, binner=\$binner"
        echo "\$binner" >> binners.txt

        # Only the samples of this batch get reports
        case " ${meta_ids} " in
            *" \$sample_id "*) ;;
            *) continue ;;
        esac

        # Linked, not copied, and named .fa whatever the extension, since
        # checkm2 runs with -x .fa
        find -L "\$bin_dir" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) | while read bin_file; do
            bin_name=\$(basename "\$bin_file")
            bin_name=\${bin_name%.*}
            ln -s "../\$bin_file" "all_bins/\${sample_id}__\${binner}__\${bin_name}.fa"
            printf '%s\\t%s\\t%s\\t%s\\n' "\${sample_id}__\${binner}__\${bin_name}" "\$sample_id" "\$binner" "\$bin_name" >> bin_map.tsv
        done
    done

    # The binner list is derived from the staged directories instead of a fixed
    # set, so no header-only reports are fabricated for binners that never ran
    # (audit #6)
    if [ ! -s binners.txt ]; then
        echo "ERROR: no binner output directories were staged/recognized" >&2
        exit 1
    fi

    total_bins=\$(( \$(wc -l < bin_map.tsv) - 1 ))
    if [ "\$total_bins" -gt 0 ]; then
        echo "Found \$total_bins genome bins across all samples and binners. Running CheckM2 assessment..."
        checkm2 predict \\
            --threads ${task.cpus} \\
            --input all_bins \\
            --output-directory checkm2_output \\
            --database_path ${checkm_db} \\
            -x .fa
        if [ ! -s checkm2_output/quality_report.tsv ]; then
            echo "ERROR: checkm2 predict wrote no checkm2_output/quality_report.tsv for \$total_bins bins" >&2
            exit 1
        fi
    else
        echo "WARNING: No genome bins found across all samples. Creating empty reports."
    fi

    # One report per sample and binner: the binner's rows of the combined
    # report, sorted by bin name, with the sample/binner prefix removed. A
    # binner without bins in any sample gets header-only reports
    python3 <<-PYEOF
	import csv

	meta_ids = "${meta_ids}".split()
	binners = sorted(set(open("binners.txt").read().split()))
	empty_header = "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes\\n"

	with open("bin_map.tsv") as f:
	    bin_map = {row["name"]: row for row in csv.DictReader(f, delimiter="\\t")}
	binners_with_bins = {row["binner"] for row in bin_map.values()}

	header, rows = None, {}
	if bin_map:
	    with open("checkm2_output/quality_report.tsv") as f:
	        header = f.readline()
	        for line in f:
	            if not line.strip():
	                continue
	            name, rest = line.split("\\t", 1)
	            entry = bin_map[name]
	            key = (entry["sample_id"], entry["binner"])
	            rows.setdefault(key, []).append((entry["bin"], entry["bin"] + "\\t" + rest))

	for binner in binners:
	    for meta_id in meta_ids:
	        with open(f"{meta_id}_{binner}_quality_report.tsv", "w") as out:
	            if binner not in binners_with_bins:
	                out.write(empty_header)
	                continue
	            out.write(header)
	            out.writelines(line for _, line in sorted(rows.get((meta_id, binner), [])))
	PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version 2>&1 | sed 's/checkm2: version //')
    END_VERSIONS
    """
    
    stub:
    def meta_ids = meta_list.collect { m -> m instanceof Map ? (m.id ?: m['id']) : m.toString() }.join(' ')
    """
    # Create stub quality report files for each sample and binner
    for meta_id in ${meta_ids}; do
        for binner in metabat semibin comebin metawrap; do
            cat > \${meta_id}_\${binner}_quality_report.tsv << 'EOF'
    Name	Completeness	Contamination
    bin.1	95.5	2.1
    bin.2	87.3	1.5
    EOF
        done
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: 1.1.0
    END_VERSIONS
    """
}
