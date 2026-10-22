/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Test harness for SOFTWARE_VERSIONS (tests/subworkflows/software_versions.nf.test)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A process that reports its versions the way current nf-core modules do
    (one tuple per tool on the 'versions' topic), so the test reads them
    through channel.topic('versions') exactly as main.nf does
----------------------------------------------------------------------------------------
*/

include { SOFTWARE_VERSIONS } from '../../subworkflows/local/software_versions'

process TOPIC_VERSIONS_EMITTER {
    tag "${id}"
    label 'process_single'
    container 'ubuntu:24.04'

    input:
    val id

    output:
    tuple val("${task.process}"), val('fastp'), eval('echo 1.4.0'), topic: versions
    tuple val("${task.process}"), val('pigz'), eval('echo 2.8'), topic: versions

    script:
    """
    true
    """
}

workflow TOPIC_VERSIONS_TEST {

    take:
    ch_ids            // channel: task ids, one task each
    ch_versions_yml   // channel: versions.yml files

    main:
    TOPIC_VERSIONS_EMITTER( ch_ids )
    SOFTWARE_VERSIONS( ch_versions_yml, channel.topic('versions') )

    emit:
    yml = SOFTWARE_VERSIONS.out.yml
}
