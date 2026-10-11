# K-MARVEL 1.0.0 for the BugBuster read-level AMR step (KMARVEL module).
#
# Upstream (https://bitbucket.org/amr-avengers/k-marvel) ships no package or
# container: it is built from source with cargo, and at run time it shells
# out to `blastp` (BLAST+ >= 2.12) on the assembled ARGs. This image pins
#   - K-MARVEL at commit 7108d2d (Cargo version 1.0.0), built with the Rust
#     toolchain upstream documents (cargo 1.93.1). Upstream's Cargo.lock at
#     that commit is stale (Cargo.toml asks for bio 3, rayon 1.12, ... and
#     `cargo build --locked` refuses it), so the build uses ./Cargo.lock: the
#     minimal re-resolution (`cargo update --workspace`, 2026-10-11). Drop it
#     when upstream commits a current lockfile
#   - BLAST+ 2.17.0, the release upstream's get_ref_db.sh installs, from
#     NCBI's versioned path (md5 published by NCBI, verified below); only
#     blastp and makeblastdb are kept
# plus procps (`ps`), which Nextflow needs for task metrics, and wget and
# bzip2 for FORMAT_KMARVEL_DB (the shared wget image cannot unpack CARD's
# .tar.bz2).
#
# The CARD reference is not baked in: FORMAT_KMARVEL_DB builds it in this
# image (bin/kmarvel_db_reformat.sh) and the KMARVEL module indexes it per
# task.
#
# Built and pushed to ghcr.io/gene2dis/bugbuster-kmarvel by
# .github/workflows/kmarvel-image.yml.

FROM rust:1.93.1-slim-bookworm@sha256:5b9332190bb3b9ece73b810cd1f1e9f06343b294ce184bcb067f0747d7d333ea AS build

ARG KMARVEL_COMMIT=7108d2d4e0116115f7b09c99ee3ac144eb25d6f4

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# make: jemalloc-sys (the jemallocator dependency) builds jemalloc from C
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends ca-certificates curl make \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /src/kmarvel \
    && curl -fsSL "https://bitbucket.org/amr-avengers/k-marvel/get/${KMARVEL_COMMIT}.tar.gz" \
        | tar -xz -C /src/kmarvel --strip-components=1
COPY Cargo.lock /src/kmarvel/Cargo.lock
RUN cd /src/kmarvel \
    && cargo build --release --locked \
    && install -m 0755 target/release/kmarvel /usr/local/bin/kmarvel

FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587

ARG KMARVEL_COMMIT=7108d2d4e0116115f7b09c99ee3ac144eb25d6f4
ARG BLAST_VERSION=2.17.0
ARG BLAST_MD5=bdec166721de3b55f90a3badc83538e8

LABEL org.opencontainers.image.title="K-MARVEL" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.revision="${KMARVEL_COMMIT}" \
      org.opencontainers.image.description="K-MARVEL 1.0.0 (commit 7108d2d) with BLAST+ 2.17.0 and procps for the BugBuster pipeline" \
      org.opencontainers.image.source="https://github.com/gene2dis/BugBuster" \
      org.opencontainers.image.licenses="GPL-3.0"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends bzip2 ca-certificates curl libgomp1 procps wget \
    && rm -rf /var/lib/apt/lists/*

# BLAST+: blastp (called by kmarvel) and makeblastdb (the KMARVEL module
# indexes the CARD FASTA per task)
RUN curl -fsSL -o /tmp/blast.tar.gz \
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" \
    && echo "${BLAST_MD5}  /tmp/blast.tar.gz" | md5sum -c - \
    && tar -xzf /tmp/blast.tar.gz -C /tmp \
        "ncbi-blast-${BLAST_VERSION}+/bin/blastp" "ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb" \
    && install -m 0755 "/tmp/ncbi-blast-${BLAST_VERSION}+/bin/blastp" "/tmp/ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb" /usr/local/bin/ \
    && rm -rf /tmp/blast.tar.gz "/tmp/ncbi-blast-${BLAST_VERSION}+"

COPY --from=build /usr/local/bin/kmarvel /usr/local/bin/kmarvel

# Smoke checks
RUN kmarvel --version \
    && blastp -version \
    && makeblastdb -version \
    && wget --version | head -1 \
    && bzip2 --help 2>&1 | head -1 \
    && ps -e -o pid= -o ppid= > /dev/null

CMD ["kmarvel", "--help"]
