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
    channel
        .fromPath(samplesheet)
        .ifEmpty { error "Cannot find samplesheet file: ${samplesheet}" }
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> validate_input(row) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
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

    // Validate sample name (no spaces, special chars)
    if (meta.id ==~ /.*\s.*/) {
        error "Invalid sample name '${meta.id}': Sample names cannot contain spaces"
    }

    // Check file existence
    def r1_file = file(row.r1.toString().trim(), checkIfExists: true)
    def r2_file = file(row.r2.toString().trim(), checkIfExists: true)

    // Validate file extensions
    def valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    if (!valid_extensions.any { r1_file.name.endsWith(it) }) {
        error "Invalid R1 file extension for sample '${meta.id}': ${r1_file.name}. Must be one of: ${valid_extensions.join(', ')}"
    }
    if (!valid_extensions.any { r2_file.name.endsWith(it) }) {
        error "Invalid R2 file extension for sample '${meta.id}': ${r2_file.name}. Must be one of: ${valid_extensions.join(', ')}"
    }

    // Handle optional singleton reads
    def reads = []
    if (row.s && row.s.toString().trim()) {
        def s_file = file(row.s.toString().trim(), checkIfExists: true)
        if (!valid_extensions.any { s_file.name.endsWith(it) }) {
            error "Invalid singleton file extension for sample '${meta.id}': ${s_file.name}. Must be one of: ${valid_extensions.join(', ')}"
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
    def duplicates = ids.groupBy { it }.findAll { it.value.size() > 1 }.keySet()
    if (duplicates) {
        error "Duplicate sample IDs found in samplesheet: ${duplicates.join(', ')}"
    }
}
