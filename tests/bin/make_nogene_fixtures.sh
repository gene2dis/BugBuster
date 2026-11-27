#!/bin/bash
#
# Regenerates the gene-less assembly of tests/data/functional/ (roadmap
# Phase 29, CL5; design doc Q24): an assembly whose contig passes
# --bbmap_length but gives Pyrodigal no gene, for
# tests/subworkflows/functional_annotation.nf.test.
#   nogene_contigs.fa  one 1,500 bp contig, TTAGTTAA repeated: a TAA or TAG
#                      stop in every frame of both strands (TGA alone is not
#                      enough: Pyrodigal's metagenomic models include
#                      translation table 4, where TGA is tryptophan), so
#                      Pyrodigal 3.7.1 (-p meta) writes an empty FAA and GFF
#   nogene.bam         its header-only BAM (no read aligned), as the counting
#                      BAM of a sample whose reads miss the contig
# Synthetic, not real data. The BAM is written with the samtools image of
# config/modules.config (--no-PG), so a rerun gives the same bytes. Needs docker.
#
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUT="$( cd "${SCRIPT_DIR}/../data/functional" && pwd )"
IMAGE="quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

seq=$(printf 'TTAGTTAA%.0s' $(seq 1 188) | cut -c1-1500)
{
    echo ">nogene_contig_1"
    echo "${seq}" | fold -w 60
} > "${OUT}/nogene_contigs.fa"

printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:nogene_contig_1\tLN:1500\n' > "${OUT}/nogene.sam"
docker run --rm -u "$(id -u):$(id -g)" -v "${OUT}:/d" -w /d "${IMAGE}" \
    samtools view --no-PG -b -o nogene.bam nogene.sam
rm "${OUT}/nogene.sam"
ls -l "${OUT}/nogene_contigs.fa" "${OUT}/nogene.bam"
