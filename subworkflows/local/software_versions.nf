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
    so the files are not indented; stripIndent() is Java's (JDK 15+), which
    counts the empty line after the final newline, so it changes nothing in
    a file ending in a newline (a broken heredoc, e.g. KARGVA's, stays as it
    is; parse_versions in bin/summary_common.py copes with it).
----------------------------------------------------------------------------------------
*/

workflow SOFTWARE_VERSIONS {

    take:
    ch_versions_yml     // channel: versions.yml files
    ch_topic_versions   // channel: [ process, tool, version ], e.g. channel.topic('versions')

    main:
    ch_topic_yml = ch_topic_versions
        .map { name, tool, version -> [ name.toString(), "    ${tool}: ${version.toString().trim()}".toString() ] }
        .unique()
        .groupTuple()
        .map { name, tools -> ( "\"${name}\":\n" + tools.sort().join('\n') + '\n' ).toString() }

    ch_yml = ch_versions_yml
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
