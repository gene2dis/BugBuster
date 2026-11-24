/*
 * Shared nf-test helpers (nf-test puts tests/lib on the test classpath).
 */
class TestData {

    // Local two-sample samplesheet: the first 500 read pairs of the nf-core
    // reads the test profile downloads (tests/bin/make_test_read_fixtures.sh).
    // GitHub raw downloads failed intermittently on CI runners (2026-10-09),
    // so no test reads them: every nf-test, including the real QC-only run,
    // takes this sheet, and the CI profile runs pass the same reads with
    // --input (roadmap Phase 25, RV27; review R8-11). Only -profile test
    // without --input still downloads assets/test_samplesheet.csv's reads.
    static String localSamplesheet(outputDir, projectDir) {
        def sheet = new File("${outputDir}", "local_test_samplesheet.csv")
        sheet.parentFile.mkdirs()
        def reads = "${projectDir}/tests/data/reads"
        sheet.text = "sample,r1,r2,s\n" +
            "test_sample1,${reads}/test_sample1_R1.fastq.gz,${reads}/test_sample1_R2.fastq.gz,\n" +
            "test_sample2,${reads}/test_sample2_R1.fastq.gz,${reads}/test_sample2_R2.fastq.gz,\n"
        return sheet.absolutePath
    }

    // Three samples (test_sample3 reuses test_sample1's reads), for the
    // co-assembly group tests: two groups, one of them a one-sample group
    static String localSamplesheet3(outputDir, projectDir) {
        def sheet = new File("${outputDir}", "local_test_samplesheet3.csv")
        sheet.parentFile.mkdirs()
        def reads = "${projectDir}/tests/data/reads"
        sheet.text = "sample,r1,r2,s\n" +
            "test_sample1,${reads}/test_sample1_R1.fastq.gz,${reads}/test_sample1_R2.fastq.gz,\n" +
            "test_sample2,${reads}/test_sample2_R1.fastq.gz,${reads}/test_sample2_R2.fastq.gz,\n" +
            "test_sample3,${reads}/test_sample1_R1.fastq.gz,${reads}/test_sample1_R2.fastq.gz,\n"
        return sheet.absolutePath
    }

    // n samples (sample1 ... sample<n>) on test_sample1's reads, for launch
    // checks that count samples (the --multisample_binning group size warning)
    static String localSamplesheetN(outputDir, projectDir, int n) {
        def sheet = new File("${outputDir}", "local_test_samplesheet_n.csv")
        sheet.parentFile.mkdirs()
        def reads = "${projectDir}/tests/data/reads"
        sheet.text = "sample,r1,r2,s\n" + (1..n).collect { i ->
            "sample${i},${reads}/test_sample1_R1.fastq.gz,${reads}/test_sample1_R2.fastq.gz,\n"
        }.join('')
        return sheet.absolutePath
    }

    // A samplesheet whose samples carry the given ids, all on test_sample1's
    // reads, for the sample-id launch checks; header is the samplesheet's
    // header line (a header with spaces around 'sample' tests the trimming)
    static String samplesheetWithIds(outputDir, projectDir, List ids, String header = "sample,r1,r2,s") {
        def sheet = new File("${outputDir}", "ids_samplesheet.csv")
        sheet.parentFile.mkdirs()
        def reads = "${projectDir}/tests/data/reads"
        sheet.text = "${header}\n" + ids.collect { id ->
            "${id},${reads}/test_sample1_R1.fastq.gz,${reads}/test_sample1_R2.fastq.gz,\n"
        }.join('')
        return sheet.absolutePath
    }

    // A --sample_metadata table with the given text (TSV, or CSV by a .csv name)
    static String metadataTable(outputDir, String name, String text) {
        def table = new File("${outputDir}", name)
        table.parentFile.mkdirs()
        table.text = text
        return table.absolutePath
    }

    // Every file a run published under outputDir, relative to it, without the
    // run-specific pipeline_info reports and the test's own samplesheets: the
    // list to compare with a tests/data/stub_paths/ file
    static List<String> publishedPaths(outputDir) {
        def out = new File("${outputDir}")
        def published = []
        out.eachFileRecurse { f ->
            if (f.isFile()) {
                def rel = out.toPath().relativize(f.toPath()).toString()
                if (!(rel ==~ /pipeline_info\/(execution_|pipeline_dag|params_).*/) && !rel.startsWith('local_test_samplesheet')) {
                    published << rel
                }
            }
        }
        return published
    }

    // The database folders a run must not publish under outputDir (roadmap
    // RV14: databases live in databases_dir only; storeDir with publishDir off).
    // Checked as folders, so an empty one (DeepARG's stub) counts too
    static List<String> databaseFolders(outputDir) {
        return ['format', 'download', 'rgi'].findAll { d -> new File("${outputDir}", d).exists() }
    }

    // The Nextflow version a pipeline test ran on, [major, minor], read from
    // the run's log (nf-test launches with -quiet, so the console has no
    // banner; the log has "Version: 26.10.5 build 14654"): for the tests
    // whose expected outcome differs by version (roadmap Phase 25, RV27:
    // nf-test runs on 25.10.7 and the newest Nextflow in CI)
    static List<Integer> nextflowVersion(outputDir) {
        def log = new File("${outputDir}/../meta/nextflow.log")
        def m = log.exists() ? (log.text =~ /(?m)^\s*Version: (\d+)\.(\d+)\.\d+ build/) : null
        if (m == null || !m.find()) {
            throw new IllegalStateException("no 'Version: <x.y.z> build' line in ${log}")
        }
        return [m.group(1) as Integer, m.group(2) as Integer]
    }

    // True when the run's Nextflow is older than major.minor
    static boolean nextflowBefore(outputDir, int major, int minor) {
        def (ma, mi) = nextflowVersion(outputDir)
        return ma < major || (ma == major && mi < minor)
    }
}
