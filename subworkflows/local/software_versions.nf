/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTWARE_VERSIONS: pipeline_info/software_versions.yml
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two sources, one file in one format:
    - versions.yml files, written by local modules (and nf-core modules still
      on the old interface);
    - [process, tool, version] tuples from channel.topic('versions'), which
      current nf-core modules emit instead of a file (one tuple per tool).
    Each tuple becomes the text a local module writes,
    "PROCESS":\n    tool: version\n, with the tools of one process in one
    block, so RESULTS_SUMMARY reads the file as before.

    Per-sample tasks of one process write byte-identical versions.yml files
    (and emit identical tuples), so content-level unique() collapses them to
    one block per process. Nextflow already strips the script's indentation,
    so the files are not indented. stripIndent() is Java's (JDK 15+): it
    trims trailing whitespace on each line, but it counts the empty line
    after the final newline, so it removes no indentation from a file ending
    in a newline (a broken heredoc, e.g. KARGVA's, stays indented;
    parse_versions in bin/summary_common.py copes with it).

    A topic tuple mixed into ch_versions by mistake (a module's
    `.out.versions_<tool>` after an nf-core update) is handled as a tuple, not
    read as a file; the topic reports it too, and unique() drops the copy. A
    multi-line version string is joined into one line to keep the file
    line-based.
----------------------------------------------------------------------------------------
*/

workflow SOFTWARE_VERSIONS {

    take:
    ch_versions_yml     // channel: versions.yml files (a stray [ process, tool, version ] tuple is accepted)
    ch_topic_versions   // channel: [ process, tool, version ], e.g. channel.topic('versions')

    main:
    ch_versions_in = ch_versions_yml
        .branch { v ->
            tuple: v instanceof List
            file: true
        }

    ch_topic_yml = ch_topic_versions
        .mix( ch_versions_in.tuple )
        .map { name, tool, version ->
            [ name.toString(), "    ${tool}: ${version.toString().trim().replaceAll(/\s*\n\s*/, ' ')}".toString() ]
        }
        .unique()
        .groupTuple()
        .map { name, tools -> ( "\"${name}\":\n" + tools.sort().join('\n') + '\n' ).toString() }

    ch_yml = ch_versions_in.file
        .map { yml -> yml.text.stripIndent() }
        .mix( ch_topic_yml )
        .mix( channel.of(
            ( "\"${workflow.manifest.name}\":\n" +
              "    pipeline: ${workflow.manifest.version}\n" +
              "    nextflow: ${nextflow.version}\n" ).toString()
        ) )
        .unique()
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.output}/pipeline_info",
            sort: true,
            newLine: false
        )

    emit:
    yml = ch_yml   // channel: software_versions.yml
}
