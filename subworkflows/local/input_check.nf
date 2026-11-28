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
    check_sample_ids(samplesheet)

    emit:
    // channel: [ val(meta), [ reads ] ] (the one emit, unnamed: INPUT_CHECK.out)
    channel
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
        error "Invalid samplesheet: 'r2' column is empty for sample: ${row.sample}. The pipeline needs paired-end reads (r1 and r2); the s column is only for singletons of a paired library"
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

    // HUMAnN 4.0.0a2 crashes on a sample name holding 's__' and 't__' (it
    // misreads the MetaPhlAn profile's header line as a taxon row, design doc
    // Q2 bug 2): such names contain '__', which check_sample_ids rejects for
    // every run (unsafeIdReason)

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
 * Get list of sample IDs from samplesheet, read as INPUT_CHECK's channel reads
 * it (fields trimmed) and with a byte order mark removed from the header, so
 * the launch checks and the group tables see the same ids as the run
 * @param samplesheet Path to samplesheet CSV
 * @return List of sample IDs
 */
def get_sample_ids(samplesheet) {
    def lines = file(samplesheet).splitCsv(sep: ',', strip: true)
    if (!lines) {
        return []
    }
    def header = lines[0].collect { c -> c.toString().trim() }
    header[0] = header[0].replaceFirst(/^\uFEFF/, '').trim()
    def col = header.indexOf('sample')
    if (col < 0) {
        return []
    }
    def ids = []
    lines.drop(1).each { fields ->
        def id = col < fields.size() ? fields[col]?.toString()?.trim() : null
        if (id) {
            ids << id
        }
    }
    return ids
}

// Ids the pipeline gives to its own pooled assemblies and groups: without
// --coassembly_group_by the co-assembly is 'coassembly', without
// --multisample_group_by the multi-sample group is 'multisample'
def reservedIds() {
    return ['coassembly', 'multisample']
}

/**
 * Why an id (a sample id, or a co-assembly / multi-sample group label, which
 * is the id of its assembly) is unsafe, or null when it is safe. Ids are read
 * back from TSV/CSV tables, where pandas and R turn a number-like value
 * (001 -> 1, 1e5 -> 100000.0, inf), a missing-value word (NA, NaN, null,
 * None) or a boolean word (true, false) into something else; GTDB_TK_BATCH
 * names pooled bins <id>__<bin> and splits them back at the first '__' (so
 * an id may not hold '__' nor end in '_'); and
 * the reserved ids name the pipeline's own pooled assemblies and groups.
 */
def unsafeIdReason(id) {
    def lower = id.toString().toLowerCase()
    if (lower ==~ /[0-9]+(\.[0-9]*)?(e[+-]?[0-9]+)?/ || lower in ['inf', 'infinity']) {
        return "table readers take it for a number"
    }
    if (lower in ['na', 'n/a', 'nan', 'null', 'none']) {
        return "table readers take it for a missing value"
    }
    if (lower in ['true', 'false']) {
        return "table readers take it for a true/false value"
    }
    if (id.toString().contains('__')) {
        return "it contains '__', which separates the sample id from the bin name in GTDB-Tk's pooled bins"
    }
    if (id.toString().endsWith('_')) {
        return "it ends in '_', so <id>__<bin> holds '___' and the id is read back without its last '_' (GTDB-Tk's pooled bins, the CheckM2 and GTDB-Tk report names)"
    }
    if (id.toString() in reservedIds()) {
        return "'${id}' is the id the pipeline gives to its own pooled ${id == 'coassembly' ? 'co-assembly' : 'multi-sample group'}"
    }
    return null
}

/**
 * Check the samplesheet's header at launch: the sample, r1 and r2 columns must
 * be there. Without this, a header such as sample_id,R1,R2 failed row by row
 * with "'sample' column is empty", which names the symptom, not the header
 * (review R3-08). An empty file is left to the "contains no samples" error
 * @param samplesheet Path to samplesheet CSV
 */
def check_samplesheet_header(samplesheet) {
    def lines = file(samplesheet).splitCsv(sep: ',', strip: true)
    if (!lines) {
        return
    }
    def header = lines[0].collect { c -> c.toString().trim() }
    header[0] = header[0].replaceFirst(/^\uFEFF/, '').trim()
    def missing = ['sample', 'r1', 'r2'].findAll { c -> !(c in header) }
    if (missing) {
        error "Invalid samplesheet header in ${samplesheet}: no ${missing.join(', ')} column${missing.size() > 1 ? 's' : ''} (header: ${header.join(',')}). The header must be sample,r1,r2,s (lower case): paired-end reads in r1 and r2, and an optional singleton FASTQ of the same library in s"
    }
}

/**
 * Check the samplesheet's sample IDs at launch: no duplicates, and no id
 * that unsafeIdReason() rejects. One error lists every offending id
 * @param samplesheet Path to samplesheet CSV
 */
def check_sample_ids(samplesheet) {
    check_samplesheet_header(samplesheet)
    def ids = get_sample_ids(samplesheet)
    def duplicates = ids.groupBy { id -> id }.findAll { entry -> entry.value.size() > 1 }.keySet()
    if (duplicates) {
        error "Duplicate sample IDs found in samplesheet: ${duplicates.join(', ')}"
    }
    def unsafe = ids.collect { id -> [id, unsafeIdReason(id)] }.findAll { pair -> pair[1] != null }
    if (unsafe) {
        error "Invalid sample name(s) in ${samplesheet}: ${unsafe.collect { id, reason -> "'${id}' (${reason})" }.join('; ')}. Rename them in the samplesheet's sample column (and in --sample_metadata) so that none reads as a number (001, 1.5, 1e5, inf), is a missing-value or true/false word (NA, NaN, null, None, true, false, in any case) or a reserved id (coassembly, multisample), or contains '__' or ends in '_': e.g. S001 for 001, site_rep1 for site__rep1, S1 for S1_"
    }
}
