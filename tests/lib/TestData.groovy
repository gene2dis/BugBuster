/*
 * Shared nf-test helpers (nf-test puts tests/lib on the test classpath).
 */
class TestData {

    // Local two-sample samplesheet for the stub and validation tests: the
    // first 500 read pairs of the nf-core reads the test profile downloads
    // (tests/bin/make_test_read_fixtures.sh). Stub runs ignore read content,
    // and GitHub raw downloads failed intermittently on CI runners
    // (2026-10-09); the real-data test still uses assets/test_samplesheet.csv.
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

    // A --sample_metadata table with the given text (TSV, or CSV by a .csv name)
    static String metadataTable(outputDir, String name, String text) {
        def table = new File("${outputDir}", name)
        table.parentFile.mkdirs()
        table.text = text
        return table.absolutePath
    }
}
