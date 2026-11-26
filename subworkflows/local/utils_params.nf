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
    return file(sourmashCustomFiles()[0]).baseName
}

// --custom_sourmash_db as its [k-mer file, lineages file] list. A params file
// or a config gives a list; on the command line the value is a String, taken
// as a JSON-style list ('["kmers.zip","lineages.csv"]') or a comma-separated
// pair ('kmers.zip,lineages.csv'). main.nf stops unless it holds two entries
// (review R1-02: the documented command-line form failed validation)
def sourmashCustomFiles() {
    def value = params.custom_sourmash_db
    if (!value) {
        return []
    }
    def entries = value instanceof List ? value : value.toString().trim().replaceAll(/^\[|\]$/, '').tokenize(',')
    return entries.collect { e -> e.toString().trim().replaceAll(/^["']|["']$/, '').trim() }.findAll { e -> e }
}

// --binners as a list of binner names, trimmed and lower-cased, each once:
// main.nf validates it and BINNING runs from it, so both read it the same way
// (roadmap Phase 20, RV4). A repeated name (metabat2,metabat2) counted as two
// binners, which turned on Binette with one bin set and passed the ">= 2
// binners" checks of --dereplicate_mags and --mag_level_functional
def parseBinners() {
    def names = params.binners instanceof List ? params.binners : params.binners.toString().tokenize(',')
    return names.collect { b -> b.toString().trim().toLowerCase() }.findAll { b -> b }.unique()
}

// true when a (plain, uncompressed) FASTA holds at least one record. BBMAP
// writes a FASTA without sequences when every contig of an assembly is below
// --bbmap_length; ASSEMBLY drops such an assembly before the contig stages
// (roadmap Phase 21, RV10), and the geNomad filter in main.nf uses it too
def hasFastaRecord(fasta) {
    return fasta.withReader { reader -> reader.readLine()?.startsWith('>') ?: false }
}

// --bbmap_length as BBMAP reads it: an Integer, 0 included (roadmap Phase 21, RV11)
def bbmapMinLength() {
    return (params.bbmap_length != null ? params.bbmap_length : 1000) as Integer
}

// The launch error for a cloud profile (aws, gcp, azure) run without an
// object-storage work directory, or null when the work directory is fine. The
// work directory may come from the profile's --*_workdir parameter or from
// Nextflow's -work-dir. A cloud Path's toString() drops the scheme
// (s3://my-bucket/work gives /my-bucket/work, gs://my-bucket/work gives
// /work), so the scheme is read from toUriString(); comparing toString() made
// -work-dir alone always fail (roadmap Phase 26, RV30)
def cloudWorkdirError(profile, workDir) {
    def requirements = [
        aws  : ['aws_workdir', 's3://my-bucket/work', 's3://'],
        gcp  : ['gcp_workdir', 'gs://my-bucket/work', 'gs://'],
        azure: ['azure_workdir', 'az://my-container/work', 'az://'],
    ]
    def active = (profile ?: '').toString().tokenize(',')*.trim()
    def uri = workDir == null ? '' : (workDir instanceof java.nio.file.Path ? workDir.toUriString() : workDir.toString())
    def failing = requirements.find { name, req -> name in active && !params[req[0]] && !uri.startsWith(req[2]) }
    if (!failing) {
        return null
    }
    def (param_name, example, _scheme) = failing.value
    return "-profile ${failing.key} requires an object-storage work directory: pass --${param_name} ${example} (or -work-dir ${example})"
}
