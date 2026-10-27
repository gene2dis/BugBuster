import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/*
 * What a task staged, read from its work directory: nf-test sees a task's
 * outputs but not its inputs, and staged inputs are links in the directory
 * that holds its outputs.
 */
class TaskInputs {

    // Names of the FASTQ files staged in the work directory of the task that
    // wrote the given output file
    static List<String> stagedReads(output) {
        return stagedReadsIn(Paths.get(output.toString()).parent)
    }

    // Names of the FASTQ files staged in the work directory of the task that
    // wrote the file a staged link in the given output's directory points to
    // (one hop upstream: e.g. BBMAP's contigs link -> the MEGAHIT task)
    static List<String> upstreamStagedReads(output, String linkName) {
        def link = Paths.get(output.toString()).parent.resolve(linkName)
        return stagedReadsIn(Files.readSymbolicLink(link).parent)
    }

    // Targets of the BAM links staged in the work directory of the task that
    // wrote the given output, in staging order (nf-core modules stage a list
    // as 1/<file>, 2/<file>, ...)
    static List<String> stagedBams(output) {
        def dir = Paths.get(output.toString()).parent
        def links = []
        Files.walk(dir, 2).withCloseable { paths ->
            paths.filter { p -> Files.isSymbolicLink(p) && p.toString().endsWith('.bam') }
                .forEach { p -> links << p }
        }
        return links
            .sort { p -> dir.relativize(p).getName(0).toString().toInteger() }
            .collect { p -> Files.readSymbolicLink(p).toString() }
    }

    static List<String> stagedReadsIn(Path dir) {
        return dir.toFile().listFiles()
            .findAll { f -> Files.isSymbolicLink(f.toPath()) && (f.name =~ /\.f(ast)?q(\.gz)?$/) }
            .collect { f -> f.name }
            .sort()
    }
}
