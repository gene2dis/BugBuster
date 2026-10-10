/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/**
 * Subworkflow for validating input samplesheet and creating read channels
 */

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    check_duplicates(samplesheet)

    emit:
    // channel: [ val(meta), [ reads ] ]
    reads = channel
        .fromPath(samplesheet)
        .ifEmpty { error "Cannot find samplesheet file: ${samplesheet}" }
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> validate_input(row) }
        // A header-only (or empty) samplesheet previously yielded a
        // "successful" run that processed nothing
        .ifEmpty { error "Samplesheet contains no samples: ${samplesheet}" }
}

/**
 * Validate input row from samplesheet
 * @param row Map with sample, r1, r2, s columns
 * @return Tuple of [meta, [reads]]
 */
def validate_input(row) {
    // Check required fields
    if (!row.sample) {
        error "Invalid samplesheet: 'sample' column is empty for row: ${row}"
    }
    if (!row.r1) {
        error "Invalid samplesheet: 'r1' column is empty for sample: ${row.sample}"
    }
    if (!row.r2) {
        error "Invalid samplesheet: 'r2' column is empty for sample: ${row.sample}"
    }

    // Create meta map
    def meta = [:]
    meta.id = row.sample.toString().trim()

    // Sample ids become file and directory names throughout the pipeline:
    // restrict to path-safe characters (previously only whitespace was
    // rejected, letting '/', ':' etc. flow into filenames)
    if (!(meta.id ==~ /[A-Za-z0-9][A-Za-z0-9_.\-]*/)) {
        error "Invalid sample name '${meta.id}': use only letters, digits, underscore, dot or hyphen, starting with a letter or digit"
    }

    // HUMAnN 4.0.0a2 parses every MetaPhlAn profile line containing both
    // 's__' and 't__' as a species row - including the header lines that carry
    // the sample name - and crashes with an IndexError (design doc Q2 bug 2,
    // fixed only in unreleased upstream code)
    if (params.read_level_functional == 'humann' && (meta.id.contains('s__') || meta.id.contains('t__'))) {
        error "Invalid sample name '${meta.id}' for --read_level_functional humann: sample names must not contain 's__' or 't__' (HUMAnN 4.0.0a2 misreads them as MetaPhlAn taxon labels and crashes); rename the sample"
    }

    // Check file existence
    def r1_file = file(row.r1.toString().trim(), checkIfExists: true)
    def r2_file = file(row.r2.toString().trim(), checkIfExists: true)

    // A copy-paste mistake (same file for both mates) would otherwise run
    // "successfully" and produce garbage
    if (r1_file == r2_file) {
        error "Invalid samplesheet: r1 and r2 are the same file for sample '${meta.id}': ${r1_file}"
    }

    // Validate file extensions
    def valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    if (!valid_extensions.any { ext -> r1_file.name.endsWith(ext) }) {
        error "Invalid R1 file extension for sample '${meta.id}': ${r1_file.name}. Must be one of: ${valid_extensions.join(', ')}"
    }
    if (!valid_extensions.any { ext -> r2_file.name.endsWith(ext) }) {
        error "Invalid R2 file extension for sample '${meta.id}': ${r2_file.name}. Must be one of: ${valid_extensions.join(', ')}"
    }

    // Handle optional singleton reads
    def reads = []
    if (row.s && row.s.toString().trim()) {
        def s_file = file(row.s.toString().trim(), checkIfExists: true)
        if (!valid_extensions.any { ext -> s_file.name.endsWith(ext) }) {
            error "Invalid singleton file extension for sample '${meta.id}': ${s_file.name}. Must be one of: ${valid_extensions.join(', ')}"
        }
        if (s_file == r1_file || s_file == r2_file) {
            error "Invalid samplesheet: singleton file duplicates r1/r2 for sample '${meta.id}': ${s_file}"
        }
        meta.single_end = false
        meta.has_singletons = true
        reads = [r1_file, r2_file, s_file]
    } else {
        meta.single_end = false
        meta.has_singletons = false
        reads = [r1_file, r2_file]
    }

    return [meta, reads]
}

/**
 * Get list of sample IDs from samplesheet
 * @param samplesheet Path to samplesheet CSV
 * @return List of sample IDs
 */
def get_sample_ids(samplesheet) {
    def ids = []
    file(samplesheet).splitCsv(header: true).each { row ->
        if (row.sample) {
            ids << row.sample.toString().trim()
        }
    }
    return ids
}

/**
 * Check for duplicate sample IDs
 * @param samplesheet Path to samplesheet CSV
 */
def check_duplicates(samplesheet) {
    def ids = get_sample_ids(samplesheet)
    def duplicates = ids.groupBy { id -> id }.findAll { entry -> entry.value.size() > 1 }.keySet()
    if (duplicates) {
        error "Duplicate sample IDs found in samplesheet: ${duplicates.join(', ')}"
    }
}
