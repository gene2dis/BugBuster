/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMPLE GROUPS FROM --sample_metadata (roadmap Phase 15, A4; Phase 18b, B4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    --coassembly_group_by <column> names a column of --sample_metadata (the
    table the run summary reads: a 'sample' column, or sample names in the
    first column; TSV, or CSV by its .csv extension, as
    bin/results_summary.py load_metadata reads it). Every distinct value is one
    co-assembly, labelled with the value's characters outside [A-Za-z0-9._-]
    replaced by '_'. Without the param every sample goes into one co-assembly,
    'coassembly'. Rows for samples that are not in the samplesheet are ignored
    (a metadata table can cover several batches).
    --multisample_group_by <column> forms the groups of multi-sample binning
    (each sample's coverage on every assembly of its group) with the same
    rules (metadataGroups); without it the run is one group, 'multisample'.
----------------------------------------------------------------------------------------
*/

include { get_sample_ids } from './input_check'

/**
 * Sample -> co-assembly label under --coassembly_group_by, [:] without it.
 * Stops the run (error) on any problem with the param, the table or the
 * labels, listing the offending samples or values.
 */
def coassemblyGroups() {
    if (!params.coassembly_group_by) {
        return [:]
    }
    if (params.assembly_mode != 'coassembly') {
        error("--coassembly_group_by ${params.coassembly_group_by} groups the samples of a co-assembly, but --assembly_mode is '${params.assembly_mode}': add --assembly_mode coassembly")
    }
    return metadataGroups('coassembly_group_by', 'coassembly', 'the ungrouped co-assembly')
}

/**
 * Sample -> group label from the --sample_metadata column that the param
 * <param_name> names (--coassembly_group_by, --multisample_group_by), [:]
 * when the param is not set. The same rules for every grouping param: a
 * missing table, column or value, or labels that clash, stop the run. A
 * label may not be a sample id or <reserved>, the label of the ungrouped
 * run (<reserved_what> in the error).
 */
def metadataGroups(param_name, reserved, reserved_what) {
    if (!params[param_name]) {
        return [:]
    }
    def column = params[param_name].toString()
    def opt = "--${param_name}"
    if (!params.sample_metadata) {
        error("${opt} ${column} names a column of --sample_metadata: pass a metadata table (TSV/CSV with a 'sample' column) as well")
    }

    // The samplesheet's sample ids, as INPUT_CHECK reads them (ids and the
    // metadata's header and ids are compared trimmed)
    def samples = get_sample_ids(file(params.input, checkIfExists: true))

    def table = file(params.sample_metadata, checkIfExists: true)
    def sep = table.name.toLowerCase().endsWith('.csv') ? ',' : '\t'
    def lines = table.splitCsv(sep: sep)
    if (!lines) {
        error("--sample_metadata ${params.sample_metadata} is empty")
    }
    def header = lines[0].collect { c -> c.toString().trim() }
    header[0] = header[0].replaceFirst(/^\uFEFF/, '').trim()
    def rows = lines.drop(1).collect { fields -> [header, fields].transpose().collectEntries { name, value -> [(name): value] } }
    def key = 'sample' in header ? 'sample' : header[0]
    if (!(column in header) || column == key) {
        error("${opt} ${column}: --sample_metadata ${params.sample_metadata} has no such column (columns: ${header.findAll { c -> c != key }.join(', ')})")
    }

    def value_of = [:]
    rows.each { row ->
        def id = row[key]?.toString()?.trim()
        if (id != null && id in samples) {
            if (value_of.containsKey(id)) {
                error("--sample_metadata ${params.sample_metadata}: sample ${id} has more than one row")
            }
            value_of[id] = row[column]?.toString()?.trim() ?: ''
        }
    }
    def no_row = samples.findAll { id -> !value_of.containsKey(id) }
    if (no_row) {
        error("${opt} ${column}: no --sample_metadata row for sample(s) ${no_row.join(', ')}")
    }
    def blank = samples.findAll { id -> value_of[id] == '' }
    if (blank) {
        error("${opt} ${column}: no value in column ${column} for sample(s) ${blank.join(', ')}")
    }

    // Labels name folders and files: path-safe, distinct, and never a sample
    // id or the ungrouped co-assembly's id
    def label_of = (value_of.values() as List).unique().collectEntries { value -> [(value): value.replaceAll(/[^A-Za-z0-9._-]/, '_')] }
    def by_label = label_of.groupBy { _value, label -> label }
    def clashes = by_label.findAll { _label, values -> values.size() > 1 }
    if (clashes) {
        error("${opt} ${column}: the values ${clashes.collect { label, values -> values.keySet().collect { v -> "'${v}'" }.join(' and ') + " both become '${label}'" }.join('; ')} (characters outside A-Z, a-z, 0-9, '.', '_' and '-' become '_'): make them differ in those characters")
    }
    // folder names that differ only in letter case are one folder on a case-insensitive filesystem
    def case_clashes = label_of.groupBy { _value, label -> label.toLowerCase() }.findAll { _key, values -> values.size() > 1 }
    if (case_clashes) {
        error("${opt} ${column}: the values ${case_clashes.collect { _key, values -> values.keySet().collect { v -> "'${v}'" }.join(' and ') }.join('; ')} give labels that differ only in letter case (${case_clashes.collect { _key, values -> values.values().collect { l -> "'${l}'" }.join(' and ') }.join('; ')}), the same folder on a case-insensitive filesystem: make them differ in more than case")
    }
    label_of.each { value, label ->
        if (!(label ==~ /[A-Za-z0-9][A-Za-z0-9_.\-]*/)) {
            error("${opt} ${column}: the value '${value}' gives the label '${label}', which does not start with a letter or digit")
        }
        if (label == reserved) {
            error("${opt} ${column}: the value '${value}' gives the label '${reserved}', the name of ${reserved_what}: rename it")
        }
        if (label in samples) {
            error("${opt} ${column}: the value '${value}' gives the label '${label}', which is also a sample id: rename it")
        }
    }
    return samples.sort(false).collectEntries { id -> [(id): label_of[value_of[id]]] }
}

/**
 * The co-assembly a sample goes into: 'coassembly' without groups, else its
 * group's label. A sample the groups do not list stops the run, rather than
 * getting a co-assembly of its own.
 */
def coassemblyLabel(groups, id) {
    if (!groups) {
        return 'coassembly'
    }
    if (!groups.containsKey(id)) {
        error("--coassembly_group_by ${params.coassembly_group_by}: sample ${id} is in no co-assembly group (groups were read for: ${groups.keySet().join(', ')})")
    }
    return groups[id]
}
