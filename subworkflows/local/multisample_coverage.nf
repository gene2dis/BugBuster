/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTI-SAMPLE COVERAGE SUBWORKFLOW (--multisample_binning, roadmap Phase 18b, B4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Per-sample assembly mode: MetaBAT2's depth table for each assembly, with
    one depth column per sample of its group (--multisample_group_by, or the
    whole run), from one extra alignment per sample.

    DESCRIPTION:
        The group's assemblies are concatenated (contigs named
        <sample>:<contig>), every sample of the group is aligned once to the
        concatenation (bwa-mem3, the filters of the per-sample alignment),
        jgi_summarize_bam_contig_depths gives one depth column per sample, and
        the table is split back per assembly. Read pairs from a genome
        assembled in several samples count on one copy only, so these depths
        are competitive: MetaBAT2's input only (docs/benchmarks/
        multisample_binning_2026.md: 0.32-0.97 of a sample's depth). The MAG
        depth (BIN_DEPTH), featureCounts and SemiBin2 keep the per-sample BAM.
        Groups of one sample are left out: their assembly is binned as without
        the option.

    INPUTS:
        reads:   [meta, reads] - the clean reads (the per-sample alignment's)
        contigs: [meta, contigs] - each sample's filtered contigs

    OUTPUTS:
        depth:    [sample id, <sample>_multisample_depth.txt]
        versions: Tool version information

    The concatenation waits for every assembly of the run (a sample QC drops
    is not known at launch, so a group's size is not); its alignments and
    depths start group by group from there.
----------------------------------------------------------------------------------------
*/

include { MULTISAMPLE_CONCAT                                } from '../../modules/local/multisample_concat/main'
include { BWAMEM3_INDEX as MULTISAMPLE_BWAMEM3_INDEX       } from '../../modules/local/bwamem3_samtools/main'
include { BWAMEM3_SAMTOOLS as MULTISAMPLE_BWAMEM3_SAMTOOLS } from '../../modules/local/bwamem3_samtools/main'
include { CALCULATE_DEPTH as MULTISAMPLE_CALCULATE_DEPTH   } from '../../modules/local/calculate_depth/main'
include { MULTISAMPLE_SPLIT_DEPTH                           } from '../../modules/local/multisample_split_depth/main'
include { multisampleSharedGroups                           } from './utils_coassembly'

workflow MULTISAMPLE_COVERAGE {
    take:
    reads   // channel: [ val(meta), [ reads ] ]
    contigs // channel: [ val(meta), path(contigs) ]

    main:
    ch_versions = channel.empty()
    // sample -> group label, groups of two or more samples only
    def shared = multisampleSharedGroups()

    // [group meta (id: label, samples: sorted ids), contigs in sample order]
    ch_group_contigs = contigs
        .filter { meta, _contigs -> shared.containsKey(meta.id) }
        .map { meta, contigs_file -> [shared[meta.id], meta.id, contigs_file] }
        .groupTuple(by: 0)
        .map { label, ids, files ->
            def members = [ids, files].transpose().sort { a, b -> a[0] <=> b[0] }
            [[id: label, samples: members.collect { m -> m[0] }], members.collect { m -> m[1] }]
        }

    MULTISAMPLE_CONCAT(ch_group_contigs)
    ch_versions = ch_versions.mix(MULTISAMPLE_CONCAT.out.versions.first())

    MULTISAMPLE_BWAMEM3_INDEX(MULTISAMPLE_CONCAT.out.fasta)
    ch_versions = ch_versions.mix(MULTISAMPLE_BWAMEM3_INDEX.out.versions.first())

    // [label, group meta, concatenation, index]
    ch_references = MULTISAMPLE_CONCAT.out.fasta
        .join(MULTISAMPLE_BWAMEM3_INDEX.out.index)
        .map { gmeta, fasta, index -> [gmeta.id, gmeta, fasta, index] }

    // Every sample of a group once to its concatenation; meta.group names it
    MULTISAMPLE_BWAMEM3_SAMTOOLS(
        reads
            .filter { meta, _files -> shared.containsKey(meta.id) }
            .map { meta, files -> [shared[meta.id], meta, files] }
            .combine(ch_references, by: 0)
            .filter { _label, meta, _files, gmeta, _fasta, _index -> meta.id in gmeta.samples }
            .map { label, meta, files, _gmeta, fasta, index -> [meta + [group: label], files, fasta, index] }
    )
    ch_versions = ch_versions.mix(MULTISAMPLE_BWAMEM3_SAMTOOLS.out.versions.first())

    // One depth table per group, one column per sample in sample order (the
    // group knows its samples, so its depths start once its BAMs are done)
    MULTISAMPLE_CALCULATE_DEPTH(
        MULTISAMPLE_BWAMEM3_SAMTOOLS.out.bam_only
            .map { meta, bam -> [meta.group, meta.id, bam] }
            .combine(ch_references.map { label, gmeta, fasta, _index -> [label, gmeta, fasta] }, by: 0)
            .map { _label, id, bam, gmeta, fasta -> [groupKey([gmeta, fasta], gmeta.samples.size()), id, bam] }
            .groupTuple()
            .map { key, ids, bams ->
                def (gmeta, fasta) = key.getGroupTarget()
                [gmeta, fasta, [ids, bams].transpose().sort { a, b -> a[0] <=> b[0] }.collect { pair -> pair[1] }]
            }
    )
    ch_versions = ch_versions.mix(MULTISAMPLE_CALCULATE_DEPTH.out.versions.first())

    MULTISAMPLE_SPLIT_DEPTH(MULTISAMPLE_CALCULATE_DEPTH.out.depth.map { gmeta, _fasta, depth -> [gmeta, depth] })
    ch_versions = ch_versions.mix(MULTISAMPLE_SPLIT_DEPTH.out.versions.first())

    emit:
    depth    = MULTISAMPLE_SPLIT_DEPTH.out.depth.flatten().map { f -> [f.name - '_multisample_depth.txt', f] } // channel: [ val(sample id), path(depth) ]
    versions = ch_versions // channel: path(versions.yml)
}
