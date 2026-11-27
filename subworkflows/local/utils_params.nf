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

// true when the summary reads argNorm's CARD mapping of ARGs-OAP or DeepARG
// hits: ARGNORM_SARG runs with ARGs-OAP, ARGNORM_DEEPARG with the contig
// DeepARG step (--contig_tax_and_arg or --arg_bin_clustering with an
// assembly). The summary then needs the CARD metadata of the K-MARVEL database
// (metadata.tsv, snps.txt) for the gene families and tiers, so PREPARE_DATABASES
// prepares it even without --kmarvel_prediction (roadmap Phase 27b, owner
// 2026-11-26)
def argNormCardNeeded() {
    def args_oap = flagOn(params.args_oap_prediction) || flagOn(params.read_arg_prediction)
    def deeparg = (flagOn(params.contig_tax_and_arg) || flagOn(params.arg_bin_clustering)) && params.assembly_mode != 'none'
    return flagOn(params.results_summary) && (args_oap || deeparg)
}

// The decontamination genome registries (config/databases.config), with the
// entries a -c file still sets under their names before 2.0.0
// (bowtie_ref_host_index, bowtie_ref_genomes_for_build) added, so such a
// config keeps working (main.nf warns at launch; roadmap Phase 29, CL2). For
// a key under both names the new name's entry is used: a config being renamed
// gets what its new block says
def hostRefDb() {
    return (params.bowtie_ref_host_index ?: [:]) + (params.host_ref_db ?: [:])
}

def phixRefDb() {
    return (params.bowtie_ref_genomes_for_build ?: [:]) + (params.phix_ref_db ?: [:])
}

// The launch error for a --custom_decontamination_index that does not hold an
// index of --read_aligner, or null when it does: a Bowtie2 index cannot serve
// bwa-mem3 (a bwa-mem2 index can), nor the reverse. One check for both
// aligners (issue #29): the index each needs, the other's index (named in a
// hint when that is what the directory holds) and the hint text
def decontaminationIndexError(aligner, index_path) {
    def bwamem3_index = { name -> name.endsWith('.bwt.2bit.64') }
    def rules = [
        bwamem3: [
            expected: bwamem3_index,
            other: { name -> name ==~ /.*\.bt2l?$/ },
            hint: " It holds a Bowtie2 index: run with --read_aligner bowtie2 to keep using it, rebuild it with bwa-mem3 (or bwa-mem2), or drop --custom_decontamination_index and pass --custom_host_fasta.",
            needs: "a bwa-mem3/bwa-mem2 index (*.bwt.2bit.64)"
        ],
        bowtie2: [
            // the reverse mix-up: the index the pipeline publishes by default
            // (bwamem3_index/) passed under the bowtie2 fallback
            expected: { name -> name ==~ /.*\.[0-9]\.bt2l?$/ },
            other: bwamem3_index,
            hint: " It holds a bwa-mem3/bwa-mem2 index: drop --read_aligner bowtie2 to use it with the default aligner, or point the option at a Bowtie2 index.",
            needs: "a Bowtie2 index (*.bt2)"
        ]
    ]
    def rule = rules[aligner.toString()]
    if ( rule == null ) {
        return "--read_aligner ${aligner} has no decontamination index check (known: ${rules.keySet().join(', ')})".toString()
    }
    def index_dir = file(index_path, checkIfExists: true)
    def index_files = index_dir.isDirectory() ? index_dir.listFiles().collect { f -> f.name } : []
    if ( index_files.any { name -> rule.expected(name) } ) {
        return null
    }
    def hint = index_files.any { name -> rule.other(name) } ? rule.hint : ""
    return "--read_aligner ${aligner} needs --custom_decontamination_index to be a directory with ${rule.needs}: ${index_path}.${hint}".toString()
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
