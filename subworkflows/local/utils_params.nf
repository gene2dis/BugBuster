/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER HELPERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    On Nextflow >= 26.04 a CLI parameter arrives as a String, even an explicit
    boolean: `--include_binning false` is the String "false", which is truthy
    in Groovy, so `if (params.include_binning)` would run the branch
    (Nextflow 25.10 still converts it to a Boolean). nf-schema validates the
    value but does not convert it, and a CLI value overrides any conversion
    written in nextflow.config, so every boolean parameter is read through
    flagOn() at the point of use (design doc Q13). Modules and config files,
    which cannot include this file, use the same expression inline:
    params.<name>.toString().toBoolean()
----------------------------------------------------------------------------------------
*/

// true for Boolean true and the String "true"; false for false, "false" and
// null (nf-schema has already rejected any other value)
def flagOn(value) {
    return value.toString().toBoolean()
}

// Label of the sylph database in output file names, the Sylph DB column of
// Reads_report.csv, plots and the run summary: the --sylph_db registry key,
// or a custom database's file name without extension (gtdb-r232-c200-dbv2),
// so a custom release is never shown as the registry default
def sylphDbLabel() {
    return params.custom_sylph_db ? file(params.custom_sylph_db).baseName : params.sylph_db
}

// Labels of a custom Kraken2 or sourmash database, as sylphDbLabel() does for
// sylph: the --kraken2_db / --sourmash_db registry key, or the custom
// database's own name, so a custom database is not shown as the registry
// default in file names, the Kraken DB / Sourmash DB columns, plots and the
// summary. Kraken2's database is a directory: its whole name (k2_standard; a
// dated name keeps its dots). sourmash's is a [k-mer file, lineages] list: the
// k-mer file's name without its extension (gtdb-reps-rs226-k31.dna)
def kraken2DbLabel() {
    return params.custom_kraken_db ? file(params.custom_kraken_db).name : params.kraken2_db
}

def sourmashDbLabel() {
    if (!params.custom_sourmash_db) {
        return params.sourmash_db
    }
    def kmer_db = params.custom_sourmash_db instanceof List ? params.custom_sourmash_db[0] : params.custom_sourmash_db
    return file(kmer_db).baseName
}
